! Use emacs f90-mode:   -*- mode:f90  -*-
PROGRAM CHK_DATA

!-------------------------------------------------------------------------------
! Simple program to dump the first few elements of the data array for each
! record of an ARL packed meteorological file. Used for diagnostic testing.
! Created: 23 Nov 1999 (RRD)
!          14 Dec 2000 (RRD) - fortran90 upgrade
!          18 Oct 2001 (RRD) - expanded grid domain
!
! $Log: chk_data.f,v $
! Revision 1.13  2013/08/01 20:57:39  trn
! Added iprint argument to unpack
!
! Revision 1.13  2013/08/01 20:57:39  trn
! Added iprint argument to unpack
!
! Revision 1.12  2013/07/31 19:45:14  trn
! First revisions for interp_receptors
!
! Revision 1.11  2013/04/01 16:44:24  trn
! Add support for large (nx,ny ge 1000) grids
!
! Revision 1.10  2013/03/12 16:35:42  trn
! Additions to support extra-long label/header for large nlvl
!
! Revision 1.9  2007/04/17 17:40:12  trn
! Added header decoding of heights
!
! Revision 1.8  2007/04/02 17:45:42  trn
! More robust error handling
!
! Revision 1.7  2006/11/16 18:59:11  trn
! Changes for decoding fcst minutes from extended header
!
! Revision 1.6  2006/10/06 14:53:49  trn
! Added option for 16-bit files, revert back to using default reals
!
! Revision 1.5  2006/10/05 18:28:38  trn
! Use real*8 and added code from pakinp.f90 to make it more similar to
! hymodelc code unpacking
!
!-------------------------------------------------------------------------------

  use unpack_subs

  implicit none

  REAL,          ALLOCATABLE :: RDATA(:,:)   
  CHARACTER(1),  ALLOCATABLE :: CPACK(:)

  CHARACTER(2)               :: cgrid
  integer                    :: knx, kny
  CHARACTER(4)               :: KVAR, MODEL 
  CHARACTER(52)              :: LABEL          
  CHARACTER(256)             :: FILE, fmt
  integer, parameter :: maxlenh=15000
  CHARACTER(maxlenh)            :: HEADER
  LOGICAL                    :: FTEST

  INTEGER  :: NEXP,i1,i2,di,j1,j2,dj, kret, k, krec, IY,IM,ID,IH,IF,KL
  integer :: IYR,IMO,IDA,IHR,IFC,nxy,len
  integer ::          ICX,       MN,                                              &
         NX,       NY,        NZ,                                              &
         K_FLAG,   LENH

  integer :: ierr, itime, icontrol, tmpctrl

  real ::  PREC,VAR1,                                                 &
         POLE_LAT, POLE_LON,  REF_LAT,                                         &
         REF_LON,  SIZE,      ORIENT,                                          &
         TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
         SYNC_LAT, SYNC_LON,  DUMMY

  integer :: iroot,ivar,nvars
  integer, allocatable :: outunits(:)
  character (len=4), allocatable :: outvars(:)
  character (len=256) :: outfname
  real :: xi_coord, xj_coord
  logical :: use_16
  integer :: kol, l, num_varb
  real :: height
  integer :: k108, recl_k108, k50
  integer :: iprint=3

  integer :: ivarout, interp_type
  integer :: xi1, xi2, xj1, xj2
  real :: xout, wt_i, wt_j

  i1=1
  di=50
  i2=-1
  j1=1
  dj=50
  j2=-1
  icontrol = 0
  use_16 = .TRUE.
! directory and file name
  WRITE(*,*)'Enter file name (absolute path):'
  READ(*,'(a)') FILE
  WRITE(*,*)'Enter T/F for whether this is a 16-bit file (RAMS, DWRF, ECX):'
  READ(*,*) use_16
  WRITE(*,'(a)')'Enter type of run:',' 0 - interactive',' 1 - time inventory only ',&
       & ' 2 - variables and levels for first time period', &
       & ' 3 - also min,max,mean,sample values, for first time period only', &
       & ' 4 - everything '
  READ(*,*) icontrol
  if (icontrol .ge. 3) then
     WRITE(*,*)'Enter i and j loop limits (first, last, increment) for printouts',&
          & '(first <0:suppress printout; last<0: go to end of array size)'
     READ(*,*) i1,i2,di,j1,j2,dj
  endif
  nvars = 0
  allocate(outvars(1),outunits(1))
  outvars(:)=''
  outunits(:)=-99
  if (icontrol .ge. 4) then
     write (*,*) 'Enter number of output variables to be interpolated'
     READ(*,*,err=10,end=10) nvars
     if (nvars .ge. 1) then
        write (*,*) 'Enter interpolation type (0-nearest neighbor; 1-bilinear):'
        READ(*,*) interp_type
        write (*,*) 'Enter xi,xj output coordinates:'
        read (*,*) xi_coord, xj_coord
        deallocate(outvars,outunits)
        allocate(outvars(nvars),outunits(nvars))
        iroot=index(file,'/',back=.true.) !zero if no '/'
        outunits(1)=20
        do ivar=1,nvars
           write (*,*) 'Enter 4-character variable name:'
           read (*,'(a)') outvars(ivar)
           if (ivar .gt. 1) outunits(ivar)=outunits(ivar-1)+1
           outfname = trim(file(iroot+1:))//'_'//outvars(ivar)
           write (*,'(a,i4,2a)') 'Writing variable no',ivar,' to ',trim(outfname)
           open(unit=outunits(ivar),file=outfname)
           write (outunits(ivar),'(a,/,2e20.10,/,7a5,a20)') &
                & trim(file),xi_coord, xj_coord, &
                & 'yr','mo','da','hr','fhr','mn','lvl',outvars(ivar)
        end do
     end if
  end if
  
10 continue

! test for meteo file existence
  INQUIRE(FILE=FILE,EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(*,*)'Unable to find file: ',trim(FILE)
     STOP
  END IF

! open file to decode the standard label (50) plus the
! fixed portion (108) of the extended header
  k108=108
  k50=50
  recl_k108 = k108+k50
  fmt='(A4,I3,I2,12F7.0,3I3,I2,I4)'
  OPEN(10,FILE=FILE,RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
  READ(10,REC=1)LABEL(1:k50),HEADER(1:k108)
  if (LABEL(11:11) .eq. 'X') then
     k50=52
     recl_k108 = k108+k50
     write (*,'(a,i5)') &
          'NOTE: Detected extra-long label with label length k50=',k50
     close(10)
     OPEN(10,FILE=FILE,RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
     READ(10,REC=1)LABEL(1:k50),HEADER(1:k108)
  endif
! decode extended portion of the header
  if (HEADER(105:105) .eq. 'X') then
     k108=110
     recl_k108 = k108+k50
     fmt='(A4,I3,I2,12F7.0,3I3,I2,1x,I5)'
     write (*,'(a,i5)') &
          'NOTE: Detected extra-long  header with standard header length k108=',k108
     close(10)
     OPEN(10,FILE=FILE,RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
     READ(10,REC=1)LABEL(1:k50),HEADER(1:k108)
  end if

! decode the standard portion of the index record
  ! From metset.f: read in cgrid
  if (k50 .eq. 50) then
     READ(LABEL,'(5I2,2X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
  else
     READ(LABEL,'(5I2,4X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
  end if
  
  WRITE(*,'(A,4I5)')'Opened file       : ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'WARNING Old format meteo data grid'
     WRITE(*,*)LABEL(1:k50)
     WRITE(*,*)HEADER(1:108)
     STOP
  END IF
     
  READ(HEADER(1:k108),fmt,iostat=ierr)                &
       MODEL,    ICX,       MN,                                              &
       POLE_LAT, POLE_LON,  REF_LAT,                                         &
       REF_LON,  SIZE,      ORIENT,                                          &
       TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
       SYNC_LAT, SYNC_LON,  DUMMY,                                           &
       NX,       NY,        NZ,                                              &
       K_FLAG,   LENH

  if (ierr .ne. 0) then
     WRITE(*,'(a)') 'ERROR: decoding first header, cannot recover'
     WRITE(*,'(2a)') 'LABEL=',LABEL(1:k50)
     WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
  end if

  ! From metset.f:
  KNX=ICHAR(CGRID(1:1))
  KNY=ICHAR(CGRID(2:2))
  IF(KNX.GE.64.OR.KNY.GE.64)THEN
     NX=(KNX-64)*1000+NX
     NY=(KNY-64)*1000+NY
!!$     GRID(KG,KT)%NUMBER=KNX*10+KNY
!!$  ELSE
!!$     READ(CGRID,'(I2)')IGRID
!!$     GRID(KG,KT)%NUMBER=IGRID
  END IF

! close file and reopen with proper length
  CLOSE (10)
  NXY = NX*NY
  if (use_16) nxy=2*nxy
  LEN = NXY+k50
  OPEN(10,FILE=FILE,RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')
  if (i1 .lt. 0) i1=nx+1
  if (j1 .lt. 0) j1=ny+1
  if (i2 .lt. 0) i2=nx
  if (j2 .lt. 0) j2=ny

! print file diagnostic
  WRITE(*,'(A,4I8)')'Grid size and lrec: ',NX,NY,NXY,LEN
  WRITE(*,'(A,I5)') 'Header record size: ',LENH

! allocate array space
  ALLOCATE (RDATA(NX,NY))   
  ALLOCATE (CPACK(NXY))

! read entire file and print headers
  KREC=1
  itime = 0
  tmpctrl=icontrol
100 READ(10,REC=KREC,ERR=800)LABEL(1:k50),(CPACK(K),K=1,NXY)
  ! From metset.f: read in cgrid
    if (k50 .eq. 50) then
       READ(LABEL(1:k50),'(6I2,a2,A4,I4,2E14.7)',iostat=ierr) IY,IM,ID,IH,IF,KL,  &
            cgrid,KVAR,NEXP,PREC,VAR1
    else
       READ(LABEL(1:k50),'(5I2,1x,i3,a2,A4,I4,2E14.7)',iostat=ierr) IY,IM,ID,IH,IF,KL,  &
            cgrid,KVAR,NEXP,PREC,VAR1
    end if
    

    if (ierr .ne. 0) then
       WRITE(*,'(a,i10,a)') 'ERROR: decoding LABEL for krec=',krec,',  going on to next record'
       WRITE(*,'(2a)') 'LABEL=',LABEL(1:k50)
       krec=krec+1
       goto 100
    end if

    if (kvar .eq. 'INDX' .or. itime .eq. 0) then
       if (lenh .gt. maxlenh .or. lenh .gt. nxy) then
          write (*,*) 'Warning: ignoring header characters from ',min(maxlenh,nxy), &
               & ' to ',lenh
          lenh=min(maxlenh,nxy)
       end if
       do k=1,lenh
          header(k:k) = cpack(k)
       enddo
       WRITE(*,'(A,1x,a)')LABEL(1:k50),header(1:k108)
! decode extended portion of the header
       READ(HEADER(1:k108),fmt,iostat=ierr)              &
            MODEL,    ICX,       MN,                                              &
            POLE_LAT, POLE_LON,  REF_LAT,                                         &
            REF_LON,  SIZE,      ORIENT,                                          &
            TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
            SYNC_LAT, SYNC_LON,  DUMMY,                                           &
            NX,       NY,        NZ,                                              &
            K_FLAG,   LENH
       if (ierr .ne. 0) then
          WRITE(*,'(a)') 'ERROR: decoding HEADER'
          WRITE(*,'(2a)') 'LABEL=',LABEL(1:k50)
          WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
       end if
       ! From metset.f:
       KNX=ICHAR(CGRID(1:1))
       KNY=ICHAR(CGRID(2:2))
       IF(KNX.GE.64.OR.KNY.GE.64)THEN
          NX=(KNX-64)*1000+NX
          NY=(KNY-64)*1000+NY
!!$       GRID(KG,KT)%NUMBER=KNX*10+KNY
!!$    ELSE
!!$       READ(CGRID,'(I2)')IGRID
!!$       GRID(KG,KT)%NUMBER=IGRID
       END IF
       if (itime .eq. 0) then
          !this code taken from metset:
          kol=k108+1
          height_loop: do l=1,nz
             READ(HEADER(KOL:KOL+7),'(F6.2,I2)',iostat=ierr) HEIGHT, NUM_VARB
             if (ierr .ne. 0) then
                write (*,*) 'Error decoding level ',l,' kol=',kol
                exit height_loop
             end if
             write (*,'(a,i4,a,g15.6,a,i4)') 'INDX: l=',l,' height=',height, &
                  & ' num_varb=',num_varb
             KOL=KOL+8
             do k=1,num_varb
                if (use_16) then
                   KOL=KOL+10
                else
                   KOL=KOL+8
                end if
             end do
          end do height_loop
       endif
       itime = itime + 1
       if (itime .gt. 1 .and. icontrol .lt. 4) tmpctrl=1
    else
       if (icontrol .eq. 0) then
          READ(*,*,END=800)
          CALL UNPACK(label,icx,mn,CPACK,RDATA,NX,NY,NEXP,PREC,VAR1,i1,i2,di,j1,j2,dj,use_16,iprint)
       else
          if (tmpctrl .ge. 2) WRITE(*,'(A)')LABEL(1:k50)
          if (tmpctrl .ge. 3) CALL UNPACK(label,icx,mn,CPACK,RDATA,NX,NY,NEXP,PREC,VAR1, &
               &    i1,i2,di,j1,j2,dj,use_16,iprint)
       endif
       if (nvars .gt. 0) then
          ivarout=0
          do ivar=1,nvars
             if (kvar .eq. outvars(ivar)) ivarout=ivar
          end do
          if (ivarout .gt. 0) then
             if (interp_type .eq. 1) then
                xi1=max(1,min(int(xi_coord),nx))
                wt_i=xi_coord - xi1
                xi2=min(xi1+1,nx)
                xj1=max(1,min(int(xj_coord),ny))
                wt_j=xj_coord - xj1
                xj2=min(xj1+1,ny)
                xout = (1.-wt_i)*(1.-wt_j)*rdata(xi1,xj1) + &
                     &     wt_i *(1.-wt_j)*rdata(xi2,xj1) + &
                     & (1.-wt_i)*    wt_j *rdata(xi1,xj2) + &
                     &     wt_i *    wt_j *rdata(xi2,xj2)
             elseif (interp_type .eq. 0) then
                xi1=max(1,min(nint(xi_coord),nx))
                xj1=max(1,min(nint(xj_coord),ny))
                xout = rdata(xi1,xj1)
             else
                write (*,*) 'invalid interp_type= ',interp_type,' Aborting'
                stop 'invalid interp_type'
             end if
             write (outunits(ivarout),'(7i5,e20.10)') IY,IM,ID,IH,IF,MN,KL,xout
          end if
       endif
    endif

    KREC=KREC+1
  GO TO 100

800 continue
  do ivar=1,nvars
     close(outunits(ivar))
  end do

  STOP

END PROGRAM chk_data

