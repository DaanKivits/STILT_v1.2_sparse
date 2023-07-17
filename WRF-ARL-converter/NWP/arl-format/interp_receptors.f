! Use emacs f90-mode:   -*- mode:f90  -*-
PROGRAM interp_receptors

  !-------------------------------------------------------------------------------
  ! Program for horizontal interpolation of 2d/3d fields from ARL file
  ! to receptors (lat/lon/time) specified in input ASCII file (one line per point)
  !
  ! Will interpolate one or more fields from file and write them out
  !
  ! Will process all time periods of the ARL file, writing out output 
  ! interpolated in time
  !
  ! $Log: interp_receptors.f,v $
  ! Revision 1.5  2013/12/18 15:42:09  mellis
  ! mods related to GEOS global grid
  !
  ! Revision 1.5  2013/12/18 15:42:09  mellis
  ! mods related to GEOS global grid
  !
  ! Revision 1.4  2013/08/12 18:30:42  trn
  ! Fixed nlev for output
  !
  ! Revision 1.3  2013/08/02 19:57:16  trn
  ! Allow for missing alt component of receptor id, output receptor ids for vertical profiles
  !
  ! Revision 1.2  2013/08/01 20:58:56  trn
  ! Numerous bug fixes and improvements
  !
  ! Revision 1.1  2013/07/31 19:45:14  trn
  ! First revisions for interp_receptors
  !
  !
  ! Initial revision: based on interp_terrain.f, Revision 1.8  2013/07/30  trn
  !
  !-------------------------------------------------------------------------------

  use interp_subs
  use unpack_subs

  implicit none

  REAL,          ALLOCATABLE :: tmp_data(:,:), RDATA2d(:,:,:),RDATAhalf(:,:,:,:),RDATAfull(:,:,:,:)
  CHARACTER(1),  ALLOCATABLE :: CPACK(:)

  CHARACTER(2)               :: cgrid
  integer                    :: knx, kny
  CHARACTER(4)               :: KVAR, MODEL 
  CHARACTER(52)              :: LABEL
  character(len=80)          :: outfmt
  real                       :: missval=-999.99
  CHARACTER(256)             :: FILE, inrecep, outrecep, fmt, cfg_name=' '
  integer, parameter :: maxlenh=15000
  CHARACTER(maxlenh)            :: HEADER
  integer :: k108, recl_k108, k50
  LOGICAL                    :: FTEST

  integer :: i, irec, ierr_interp
  integer :: iu_dat=10, iu_recin=11, iu_recout=12, iu_cfg=13
  integer :: file_times(6,3), file_diffs(2)

  INTEGER  :: NEXP,i1,i2,di,j1,j2,dj, kret, k, krec, IY,IM,ID,IH,IF,KL
  integer :: IYR,IMO,IDA,IHR,IFC,nxy,len, nlev_half
  integer ::          ICX,       MN,                                              &
       NX,       NY,        NZ,                                              &
       K_FLAG,   LENH

  integer :: ierr, itime, ntimes, max_time, itime_prev, itim_indx
  integer :: itime_needed, itime_data, nmissing, nbadtim, nbadll

  real ::  PREC,VAR1,                                                 &
       POLE_LAT, POLE_LON,  REF_LAT,                                         &
       REF_LON,  SIZE,      ORIENT,                                          &
       TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
       SYNC_LAT, SYNC_LON,  DUMMY

  real, allocatable :: sigma(:)
  integer :: ivar,nvar2d,nvarhalf,nvarfull,scan,ivar2d,ivarhalf,ivarfull
  character (len=4), allocatable :: var2d(:), varhalf(:), varfull(:)
  character (len=256) :: outfname
  character (len=80), allocatable :: receptors(:)
  real, allocatable :: rec_lon(:), rec_lat(:)
  integer, allocatable :: rec_tim(:), rec_ind(:), rec_idate(:,:), rec_out(:), rec_alt(:)
  real :: xi_coord, xj_coord
  logical :: use_16, continue_reading
  logical, allocatable :: found2d(:), foundhalf(:), foundfull(:)
  integer :: kol, l, num_varb
  real :: height, wgt

  integer :: nreceptors,nrec_out
  real, allocatable :: rec_dat2d(:,:), rec_dathalf(:,:,:), rec_datfull(:,:,:)
  real, allocatable :: rec2_dat2d(:,:), rec2_dathalf(:,:,:), rec2_datfull(:,:,:)
  integer :: iarg, narg, iargc, iprint, iline, interp_type
  logical                    :: keep_missing=.FALSE.
  CHARACTER(256)             :: option, fname
  logical :: lhelp
  real :: bufxkm, bufykm, bufix, bufjy

  i1=1
  di=50
  i2=-1
  j1=1
  dj=50
  j2=-1
  use_16 = .FALSE.

  !=>check for command line arguments
  NARG=IARGC()

  lhelp = narg .eq. 0

  !     process options
  iarg = 0
  FILE=' '
  inrecep=' '
  outrecep=' '
  iprint=2
  bufxkm=0.0     ! For backwards compatibility, use no buffer by default
  bufykm=0.0
  interp_type=0  ! default to nearest neighbor

  do while (iarg .lt. narg .and. .not. lhelp)
     iarg=iarg+1
     CALL GETARG(iarg,option)
     if (option(1:1) .ne. '-') exit
     iarg=iarg+1
     CALL GETARG(iarg,fname)
     SELECT CASE (option)
     CASE ('-A','-a')
        FILE = fname
     CASE ('-B','-b')
        read(fname,*,iostat=ierr) interp_type
        if (ierr .ne. 0) then
           write (*,*) 'Error decoding interp_type from -B option:',trim(fname)
           lhelp = .TRUE.
        endif
     CASE ('-C','-c')
        cfg_name = fname
     CASE ('-I','-i')
        inrecep = fname
     CASE ('-K','-k')
        read(fname,*,iostat=ierr) missval
        keep_missing = .TRUE.
     CASE ('-O','-o')
        outrecep = fname
     CASE ('-P','-p')
        SELECT CASE (fname)
        CASE ('16')
           use_16 = .TRUE.
        CASE ('8')
           use_16 = .FALSE.
        CASE DEFAULT
           write (*,*) 'Invalid value for -P option:',trim(fname)
           lhelp=.TRUE.
        end SELECT
     CASE ('-V','-v')
        read(fname,*,iostat=ierr) iprint
        if (ierr .ne. 0) then
           write (*,*) 'Error decoding iprint from -V option:',trim(fname)
           lhelp = .TRUE.
        endif
     CASE ('-X','-x')
        read(fname,*,iostat=ierr) bufxkm
     CASE ('-Y','-y')
        read(fname,*,iostat=ierr) bufykm
     CASE DEFAULT
        write (*,*) 'Invalid option:',trim(option)
        lhelp=.TRUE.
     END SELECT
  end do

  lhelp = lhelp .or. FILE .eq. ' ' .or. inrecep .eq. ' ' .or. cfg_name .eq. ' ' &
       .or. outrecep .eq. ' ' 
  lhelp = lhelp .or. iarg .gt. narg
  write (*,'(a)') 'Running interp_receptors $Revision: 1.5 $ with the following options:'
  write (*,'(a,t50,a)') 'REQUIRED: Input ARL file [-A]:',trim(FILE)
  write (*,'(a,t50,a)') 'REQUIRED: Input configuration file [-C]:',trim(cfg_name)
  write (*,'(a,t50,a)') 'REQUIRED: Input receptor file [-I]:',trim(inrecep)
  write (*,'(a,t50,a)') 'REQUIRED: Output receptor file [-O]:',trim(outrecep)
  write (*,'(a,t50,l1)') 'Flag for 16-bit ARL precision [-P 16/-P 8]: ',use_16
  write (*,'(a,t50,i10)') 'Printout control [-V]:',iprint
  write (*,'(a,t50,i10)') 'Interpolation type [-B, 0:nearest neighbor, 1:bilinear]:',interp_type
  if (keep_missing) then
     write (*,'(a,t50,g15.6)') 'Missing value indicator [-K]:',missval
  else
     write (*,'(a)') 'Not keeping receptors outside domain/time window, use -K missval to change this'
  end if
  write (*,'(a,t50,g15.6)') 'Buffer zone (km) in x-direction [-X]:',bufxkm
  write (*,'(a,t50,g15.6)') 'Buffer zone (km) in y-direction [-Y]:',bufykm

  if (iprint .lt. 2) i1 = -1

  if (lhelp) then
     write (*,'(a)') 'Usage: '
     write (*,'(2a)') '  interp_receptors -A arl_filename -I in_receptors_filename -C config_filename', &
          &' -O out_receptors_filename [-P 16/8] [-V iprint] [-B 0/1] [-K missval] [-X bufxkm] [-Y bufykm]'
     stop 'Usage'
  endif

!DEBUG: WITH THIS SECTION OF CODE HERE, ABORTS ON CLOSE FOR iu_cfg
  ! open and read input receptor file
  if (iprint .gt. 0) write (*,*) 'Reading input receptors from file: ',trim(inrecep)
  open(unit=iu_recin,file=trim(inrecep),status='old',iostat=ierr)
  if (ierr .ne. 0) then
     write (*,*) 'Error opening inrecep file: ',trim(inrecep),' ierr=',ierr
     Stop 'error opening inrecep file'
  end if
  read(iu_recin,*) nreceptors
  allocate(receptors(nreceptors),rec_out(nreceptors), &
       rec_lon(nreceptors),rec_lat(nreceptors), rec_alt(nreceptors), &
       rec_tim(nreceptors),rec_ind(nreceptors),rec_idate(5,nreceptors))
!!$  if (nreceptors .gt. max_receptors) stop 'nreceptors too large'
!!$  allocate(rec_lon(nreceptors),rec_lat(nreceptors), &
!!$       rec_tim(nreceptors),rec_ind(nreceptors),rec_idate(5,nreceptors),rec_alt(nreceptors))
  do i=1,nreceptors
     read(iu_recin,'(a)') receptors(i)
  end do
  close(iu_recin)
  if (iprint .ge. 3) write (*,'(a,i10,a/,(a))') 'Read in ',nreceptors,' receptors:', &
       (receptors(i),i=1,nreceptors)
  call id2pos(receptors,nreceptors,rec_lon,rec_lat,rec_alt,rec_idate,ierr)
  if (ierr .ne. 0) then
     write (*,*) 'id2pos returned with error code=',ierr
     if (ierr .gt. 0) then
        stop 'Fatal id2pos error, refer to comments in unpack_subs.f for error codes'
        stop 'Fatal id2pos error'
     else
        if (ierr .eq. -1) write (*,*) ' id2pos warning: extraneous delimiters, ignored'
        if (ierr .eq. -2) write (*,*) ' id2pos warning: no alt field in receptor ids, ok (not used here)'
     end if
  end if
  
  ! open output receptor file
  if (iprint .gt. 0) write (*,*) 'Writing output receptors to file: ',trim(outrecep)
  open(unit=iu_recout,file=trim(outrecep),status='unknown',iostat=ierr)
  if (ierr .ne. 0) then
     write (*,*) 'Error opening outrecep file: ',trim(outrecep),' ierr=',ierr
     Stop 'error opening outrecep file'
  end if

  ! open and read input configuration file
  if (iprint .gt. 0) write (*,*) 'Reading input configuration (variables) from file: ',trim(cfg_name)
  open(unit=iu_cfg,file=trim(cfg_name),status='old',iostat=ierr)
  if (ierr .ne. 0) then
     write (*,*) 'Error opening configuration file: ',trim(cfg_name),' ierr=',ierr
     Stop 'error opening configuration file'
  end if
  read(iu_cfg,*) nvar2d,nvarhalf,nvarfull
  if (iprint .ge. 1) write (*,'(a,3i10)') 'Specified nvar2d,nvarhalf,nvarfull=', &
          nvar2d,nvarhalf,nvarfull
  if (nvar2d .gt. 0) then
!!$     if (nvar2d .gt. max_allvars) stop 'nvar2d too large'
     allocate(var2d(nvar2d),found2d(nvar2d))
     found2d(:)=.FALSE.
     do i=1,nvar2d
        read(iu_cfg,'(a)') var2d(i)
     end do
     if (iprint .ge. 1) write (*,'(a/,(a))') 'var2d=',(var2d(i),i=1,nvar2d)
  end if
  if (nvarhalf .gt. 0) then
!!$     if (nvarhalf .gt. max_allvars) stop 'nvarhalf too large'
     allocate(varhalf(nvarhalf),foundhalf(nvarhalf))
     foundhalf(:)=.FALSE.
     do i=1,nvarhalf
        read(iu_cfg,'(a)') varhalf(i)
     end do
     if (iprint .ge. 1) write (*,'(a/,(a))') 'varhalf=',(varhalf(i),i=1,nvarhalf)
  end if
  if (nvarfull .gt. 0) then
!!$     if (nvarfull .gt. max_allvars) stop 'nvarfull too large'
     allocate(varfull(nvarfull),foundfull(nvarfull))
     foundfull(:)=.FALSE.
     do i=1,nvarfull
        read(iu_cfg,'(a)') varfull(i)
     end do
     if (iprint .ge. 1) write (*,'(a/,(a))') 'varfull=',(varfull(i),i=1,nvarfull)
  end if
  if (iprint .ge. 1) write (*,'(a,3i10)') 'After var read nvar2d,nvarhalf,nvarfull=', &
          nvar2d,nvarhalf,nvarfull
! DEBUG: PROGRAM ABORTS ON THIS CALL IF THIS SECTION FOLLOWS THE RECEPTOR SECTION
  close(iu_cfg,iostat=ierr)
  if (ierr .ne. 0) then
     write (*,*) 'Error closing configuration file: ',trim(cfg_name),' ierr=',ierr
     Stop 'error closing configuration file'
  end if

  ! test for meteo file existence
  INQUIRE(FILE=trim(FILE),EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(*,*)'Unable to find file: ',trim(FILE)
     STOP 'Missing ARL file'
  END IF

  ! open file to decode the standard label (50) plus the
  ! fixed portion (108) of the extended header
  k108=108
  k50=50
  recl_k108 = k108+k50
  fmt='(A4,I3,I2,12F7.0,3I3,I2,I4)'
  OPEN(iu_dat,FILE=trim(FILE),RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
  READ(iu_dat,REC=1)LABEL(1:k50),HEADER(1:k108)
  if (LABEL(11:11) .eq. 'X') then
     k50=52
     recl_k108 = k108+k50
     if (iprint .gt. 0) write (*,'(a,i5)') &
          'NOTE: Detected extra-long label with label length k50=',k50
     close(iu_dat)
     OPEN(iu_dat,FILE=trim(FILE),RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
     READ(iu_dat,REC=1)LABEL(1:k50),HEADER(1:k108)
  endif
! decode extended portion of the header
  if (header(105:105) .eq. 'X') then
     k108=110
     fmt='(A4,I3,I2,12F7.0,3I3,I2,1x,I5)'
     recl_k108 = k108+k50
     if (iprint .gt. 0) write (*,'(a,i5)') &
          'NOTE: Detected extra-long header with standard header length k108=',k108
     close(iu_dat)
     OPEN(iu_dat,FILE=trim(FILE),RECL=recl_k108,ACCESS='DIRECT',FORM='UNFORMATTED')
     READ(iu_dat,REC=1)LABEL(1:K50),HEADER(1:k108)
  end if

  ! decode the standard portion of the index record
  ! From metset.f: read in cgrid
  if (k50 .eq. 50) then
     READ(LABEL,'(5I2,2X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
  else
     READ(LABEL,'(5I2,4X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
  end if
  if (iprint .gt. 0) WRITE(*,'(A,4I5)')'Opened file       : ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'ERROR: Old format meteo data grid'
     WRITE(*,*)LABEL(1:K50)
     WRITE(*,*)HEADER(1:k108)
     STOP 'Old format meteo data grid'
  END IF

  ! decode extended portion of the header
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
     WRITE(*,'(2a)') 'LABEL=',LABEL(1:K50)
     WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
     STOP 'ERROR: decoding first header, cannot recover'
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
  CLOSE (iu_dat)
  NXY = NX*NY
  if (use_16) nxy=2*nxy
  LEN = NXY+k50
  OPEN(iu_dat,FILE=trim(FILE),RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')
  if (i1 .lt. 0) i1=nx+1
  if (j1 .lt. 0) j1=ny+1
  if (i2 .lt. 0) i2=nx
  if (j2 .lt. 0) j2=ny

  ! print file diagnostic
  nlev_half=nz-1
  if (iprint .gt. 0) WRITE(*,'(A,6I8)')'Grid size and lrec, nz, nlev_half: ',NX,NY,NXY,LEN,nz,nlev_half
  if (iprint .gt. 0) WRITE(*,'(A,I5)') 'Header record size: ',LENH

  ! allocate array space
  ALLOCATE (CPACK(NXY))
  if (iprint .gt. 2) write (*,*) 'Allocated CPACK array'
  tmp_data(:,:)=0. ! needed for dummy interp_2d call in scan=1 below
  if (iprint .ge. 1) write (*,'(a,3i10)') 'Before rec_dat allocate nvar2d,nvarhalf,nvarfull=', &
          nvar2d,nvarhalf,nvarfull
  if (nvar2d .gt. 0) then
     allocate(rec_dat2d(nvar2d,nreceptors),rec2_dat2d(nvar2d,nreceptors))
     rec_dat2d=missval
  end if
  if (nvarhalf .gt. 0) then
     allocate(rec_dathalf(nlev_half,nvarhalf,nreceptors),rec2_dathalf(nlev_half,nvarhalf,nreceptors))
     rec_dathalf=missval
  end if
  if (nvarfull .gt. 0) then
     allocate(rec_datfull(0:nlev_half,nvarfull,nreceptors),rec2_datfull(0:nlev_half,nvarfull,nreceptors))
     rec_datfull=missval
  end if
  if (iprint .gt. 2) write (*,*) 'Allocated rec_dat arrays'
  ALLOCATE (tmp_data(nx,ny))
  if (iprint .gt. 2) write (*,*) 'Allocated tmp_data array'
  if (nvar2d .gt. 0) then
     ALLOCATE (RDATA2d(NX,NY,nvar2d))
     if (iprint .gt. 2) write (*,*) 'Allocated RDATA2d arrays'
  end if
  if (nvarhalf .gt. 0) then
     ALLOCATE (RDATAhalf(NX,NY,NLEV_HALF,nvarhalf))
     if (iprint .gt. 2) write (*,*) 'Allocated RDATAhalf arrays'
  end if
  if (nvarfull .gt. 0) then
     ALLOCATE (RDATAfull(NX,NY,0:NLEV_HALF,nvarfull))
     if (iprint .gt. 2) write (*,*) 'Allocated RDATAfull arrays'
  end if
  

  ! read entire file:
  ! first time through: find all times and variables
  ! find time indices for receptors
  SCAN_LOOP: do scan=1,2
     KREC=1
     itime = 0
     itime_prev = 0
     continue_reading = .TRUE.
     if (iprint .gt. 0) write (*,*) 'Starting scan= ',scan
     READ_LOOP: do while (continue_reading)
        READ(iu_dat,REC=KREC,iostat=ierr)LABEL(1:K50),(CPACK(K),K=1,NXY)
        if (ierr .ne. 0) then
           if (iprint .gt. 0) write (*,*) 'Read error ierr= ',ierr,' at krec= ',krec
           continue_reading=.FALSE.
           cycle READ_LOOP
        end if

        if (k50 .eq. 50) then
           READ(LABEL(1:k50),'(6I2,a2,A4,I4,2E14.7)',iostat=ierr) IY,IM,ID,IH,IF,KL,  &
                cgrid,KVAR,NEXP,PREC,VAR1
        else
           READ(LABEL(1:k50),'(5I2,1x,i3,a2,A4,I4,2E14.7)',iostat=ierr) IY,IM,ID,IH,IF,KL,  &
                cgrid,KVAR,NEXP,PREC,VAR1
        end if

        if (ierr .ne. 0) then
           if (iprint .gt. 0) then
              WRITE(*,'(a)') 'ERROR: decoding LABEL, going on to next record'
              WRITE(*,'(2a)') 'LABEL=',LABEL(1:K50)
              WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
           end if
           krec=krec+1
           cycle READ_LOOP
        end if

        if (kvar .eq. 'INDX' .or. itime .eq. 0) then
           ! decode time etc
           if (lenh .gt. maxlenh .or. lenh .gt. nxy) then
              if (iprint .gt. 0 .and. scan .eq. 1) &
                   write (*,*) 'Warning: ignoring header characters from ', &
                   min(maxlenh,nxy), ' to ',lenh
              lenh=min(maxlenh,nxy)
           end if
           do k=1,lenh
              header(k:k) = cpack(k)
           enddo
           if (scan .eq. 1 .and. ((itime .eq. 0 .and. iprint .gt. 0) .or. &
                iprint .gt. 1)) WRITE(*,'(A,1x,a)') LABEL(1:K50),header(1:k108)
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
              if (iprint .gt. 0) then
                 WRITE(*,'(a)') 'ERROR: decoding HEADER, going on to next record'
                 WRITE(*,'(2a)') 'LABEL=',LABEL(1:K50)
                 WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
              end if
              krec=krec+1
              cycle READ_LOOP
           end if
           ! From metset.f:
           KNX=ICHAR(CGRID(1:1))
           KNY=ICHAR(CGRID(2:2))
           IF(KNX.GE.64.OR.KNY.GE.64)THEN
              NX=(KNX-64)*1000+NX
              NY=(KNY-64)*1000+NY
!!$           GRID(KG,KT)%NUMBER=KNX*10+KNY
!!$        ELSE
!!$           READ(CGRID,'(I2)')IGRID
!!$           GRID(KG,KT)%NUMBER=IGRID
           END IF
           itime = itime + 1
           ! if first scan: store level and time info
           if (scan .eq. 1) then
              if (itime .eq. 1) then
                 ! compute buffer zone:
                 if(size .ne. 0) then
                    bufix=bufxkm/size
                    bufjy=bufykm/size
                    if (iprint .gt. 0) write(*,*) 'Based on grid spacing of size=',size,', using:', &
                         ' bufix=',bufix,', bufjy=',bufjy
                 else
                    if (iprint .eq. 0) write(*,*) 'LATLON grid (size=0), no buffer zone defined.'
                    bufix=0.0
                    bufjy=0.0
                 end if
                 allocate(sigma(nz))
                 !this code taken from metset:
                 kol=k108+1
                 height_loop: do l=1,nz
                    READ(HEADER(KOL:KOL+7),'(F6.2,I2)',iostat=ierr) sigma(l), NUM_VARB
                    if (ierr .ne. 0) then
                       if (iprint .gt. 0) write (*,*) 'Error decoding level ',l,' kol=',kol
                       exit height_loop
                    end if
                    if (iprint .ge. 2) write (*,'(a,i4,a,g15.6,a,i4)') 'INDX: l=',l,' height=',sigma(l), &
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
                 write (iu_recout,'(a,i8)') 'nsigma=',nz
                 write (iu_recout,'(g15.6)') (sigma(l),l=1,nz)
              end if
              itim_indx=min(itime,3)
              file_times(1,itim_indx)=iy
              file_times(2,itim_indx)=im
              file_times(3,itim_indx)=id
              file_times(4,itim_indx)=ih
              file_times(5,itim_indx)=mn
              file_times(6,itim_indx)=itime
           end if
        else !not an index variable
           ! check for variables
           ivar2d=0
           ivarhalf=0
           ivarfull=0
           if (iprint .gt. 2) WRITE(*,'(A)')LABEL(1:K50)
           do ivar=1,nvar2d
              if (kvar .eq. var2d(ivar)) ivar2d=ivar
           end do
           do ivar=1,nvarhalf
              if (kvar .eq. varhalf(ivar)) ivarhalf=ivar
           end do
           do ivar=1,nvarfull
              if (kvar .eq. varfull(ivar)) ivarfull=ivar
           end do
           if (scan .eq. 1) then
              if (ivar2d .gt. 0) found2d(ivar2d)=.TRUE.
              if (ivarhalf .gt. 0) foundhalf(ivarhalf)=.TRUE.
              if (ivarfull .gt. 0) foundfull(ivarfull)=.TRUE.
           else
              if (itime .ne. itime_prev) then
                 if (itime .gt. 1 .and. itime_data .ne. 0) then
                    ! Process the unpacked data from the last processed time period:
                    ! Interpolate data valid at itime_data to receptor locations
                    call interp_data
                    ! reset itime_data
                    itime_data=0
                 end if
                 ! update itime_prev, itime_needed
                 itime_prev = itime
                 continue_reading=itime .le. max_time
                 itime_needed=any(rec_ind(rec_out(1:nrec_out)) .eq. itime) .or. &
                      any(rec_ind(rec_out(1:nrec_out)) .eq. itime-1)
              endif
              ! if this time level needed for any receptor:
              if (itime_needed) then
                 itime_data=itime
                 if (ivar2d .gt. 0) &
                      CALL UNPACK(label,icx,mn,CPACK,RDATA2d(1,1,ivar2d),NX,NY,NEXP,PREC,VAR1, &
                      &    i1,i2,di,j1,j2,dj, &
                      &    use_16,iprint)
                 if (ivarhalf .gt. 0) &
                      CALL UNPACK(label,icx,mn,CPACK,RDATAhalf(1,1,kl,ivarhalf),NX,NY,NEXP,PREC,VAR1, &
                      &    i1,i2,di,j1,j2,dj, &
                      &    use_16,iprint)
                 if (ivarfull .gt. 0) &
                      CALL UNPACK(label,icx,mn,CPACK,RDATAfull(1,1,kl,ivarfull),NX,NY,NEXP,PREC,VAR1, &
                      &    i1,i2,di,j1,j2,dj, &
                      &    use_16,iprint)
              end if
           end if
        end if
        KREC=KREC+1
     end do READ_LOOP
     if (scan .eq. 1) then
        ntimes=itime
        nmissing = 0
        if (nvar2d .gt. 0) nmissing=nmissing+count(.not.found2d)
        if (nvarhalf .gt. 0) nmissing=nmissing+count(.not.foundhalf)
        if (nvarfull .gt. 0) nmissing=nmissing+count(.not.foundfull)
        if (nmissing .gt. 0) then
           if (iprint .gt. 0) &
                write (*,*) 'Some of the requested variables in ',trim(cfg_name),' not found in input'
           if (keep_missing) then
              if (iprint .gt. 0) write (*,*) 'writing out missing value indicator:',missval
           else
              write (*,*) 'Aborting'
              stop 'Some of the requested variables not found in input, aborting'
           end if

        end if
        !find time indices for receptors
        call idate_timdiff(file_times(1:5,1),file_times(1:5,2:3),file_diffs)
        call idate_timdiff(file_times(1:5,1),rec_idate,rec_tim)
        if (iprint .gt. 1) then
           write (*,'(a,i10,a,/,(6a10),/,(6i10))') 'Read ',ntimes,' time periods:', &
                'year','month','day','hour','minute','index',file_times(:,:)
           write (*,*) 'file_diffs (seconds):',file_diffs
        end if
        
        rec_ind(:)=1+rec_tim(:)/file_diffs(1) !time index for time level 1
        
        where (rec_ind .eq. ntimes) rec_ind=ntimes-1
        ! find receptors for which to interpolate data
        nrec_out=0
        nbadtim=0
        nbadll=0
        rec_out(:)=0
        max_time=0
        do i=1,nreceptors
           ! check for time in file
           if (rec_ind(i) .ge. 1 .and. rec_ind(i) .le. ntimes-1) then
              ! check for point in domain
              call interp_2d(rec_lon(i),rec_lat(i),tmp_data,nx,ny,header,maxlenh, &
                   dummy,ierr_interp,iprint,bufix,bufjy,interp_type,label)
              if (ierr_interp .eq. 0) then
                 nrec_out=nrec_out+1
                 rec_out(nrec_out)=i
                 max_time=max(max_time,rec_ind(i)+1)
              else
                 nbadll=nbadll+1
              end if
           else
              nbadtim=nbadtim+1
           end if
        end do
        if (iprint .ge. 1) write (*,'(a,i8,a,i8,a)') &
             ' Retaining ',nrec_out,' of the total ',nreceptors,' receptors', &
             ' Removed ',nbadtim,' not in time window, and ',nbadll,' outside domain.', &
             ' Max_time for scan 2= ',max_time
     else
        ! scan 2: interpolate the last time in file to receptor locations if needed
        if (itime_data .ne. 0) then
           ! Interpolate data valid at itime_data to receptor locations
           call interp_data
        end if
     endif
     if (nrec_out .lt. 1) exit SCAN_LOOP
  end do SCAN_LOOP

  close(iu_dat)

  if (keep_missing) then
     nrec_out=nreceptors
     do i=1,nreceptors
        rec_out(i)=i
     end do
  end if
  write (iu_recout,'(a,i10)') 'nvar2d=',nvar2d
  write (iu_recout,'(a,i10)') 'nvarhalf=',nvarhalf
  write (iu_recout,'(a,i10)') 'nvarfull=',nvarfull
  write (iu_recout,'(a,i10)') 'nreceptors=',nrec_out
  if (nrec_out .gt. 0) then
     ! do time interpolation

     do i=1,nrec_out
        irec=rec_out(i)
        wgt=1.+ real(rec_tim(irec))/real(file_diffs(1))-real(rec_ind(irec))
        if (nvar2d .gt. 0) then
           do ivar=1,nvar2d
              rec_dat2d(ivar,irec)=(1.-wgt)*rec_dat2d(ivar,irec) + wgt*rec2_dat2d(ivar,irec)
           end do
        end if
        if (nvarhalf .gt. 0) then
           do ivar=1,nvarhalf
              do k=1,nlev_half
                 rec_dathalf(k,ivar,irec)=(1.-wgt)*rec_dathalf(k,ivar,irec) + wgt*rec2_dathalf(k,ivar,irec)
              end do
           end do
        end if
        if (nvarfull .gt. 0) then
           do ivar=1,nvarfull
              do k=0,nlev_half
                 rec_datfull(k,ivar,irec)=(1.-wgt)*rec_datfull(k,ivar,irec) + wgt*rec2_datfull(k,ivar,irec)
              end do
           end do
        end if
     end do

     ! write out data: 2d-table of 2d data for all receptors first
     write (outfmt,'(a,i2,a,i2,a)') '(a',len_trim(receptors(1)),',',nvar2d+3,'a15)'
     write (iu_recout,outfmt) 'Recep','Lat','Lon',(var2d(i),i=1,nvar2d)
     write (outfmt,'(a,i2,a)') '(a,',nvar2d+2,'g15.6)'
     do i=1,nrec_out
        irec=rec_out(i)
        write (iu_recout,outfmt) trim(receptors(irec)),rec_lat(irec),rec_lon(irec),&
             (rec_dat2d(ivar,irec),ivar=1,nvar2d)
     end do

     do i=1,nrec_out
        irec=rec_out(i)
        ! then, for each receptor, 2d-tables of vertical profiles
        write (iu_recout,'(2a)') 'Tables of vertical profiles for: ',trim(receptors(irec))
        if (nvarhalf .gt. 0) then
           write (outfmt,'(a,i2,a)') '(',nvarhalf,'a15)'
           write (iu_recout,outfmt) (varhalf(k),k=1,nvarhalf)
           write (outfmt,'(a,i2,a)') '(',nvarhalf,'g15.6)'
           do k=1,nlev_half
              write (iu_recout,outfmt) (rec_dathalf(k,ivar,irec),ivar=1,nvarhalf)
           end do
        end if

        if (nvarfull .gt. 0) then
           write (outfmt,'(a,i2,a)') '(',nvarfull,'a15)'
           write (iu_recout,outfmt) (varfull(k),k=1,nvarfull)
           write (outfmt,'(a,i2,a)') '(',nvarfull,'g15.6)'
           do k=0,nlev_half
              write (iu_recout,outfmt) (rec_datfull(k,ivar,irec),ivar=1,nvarfull)
           end do
        end if
     end do
  end if !endif nrec_out .gt. 0

  close(iu_recout)

  write (*,*) 'Normal Program Stop'

contains 

  subroutine interp_data

    ! internal subroutine, requires no arguments

    ! interpolate data from time level itime_data to receptor locations and 
    ! store in the appropriate receptor arrays

    integer :: i, irec, k, nerr !local variables

    nerr=0
    do i=1,nrec_out
       irec=rec_out(i)
       if (rec_ind(irec) .eq. itime_data) then
          ! time level 1
          if (nvar2d .gt. 0) then
             do ivar=1,nvar2d
                call interp_2d(rec_lon(irec),rec_lat(irec),rdata2d(1,1,ivar),nx,ny, &
                     header,maxlenh, rec_dat2d(ivar,irec),ierr_interp,iprint,bufix, &
                     bufjy,interp_type,label)
                if (ierr_interp .ne. 0) nerr=nerr+1
             end do
          end if
          if (nvarhalf .gt. 0) then
             do ivar=1,nvarhalf
                do k=1,nlev_half
                   call interp_2d(rec_lon(irec),rec_lat(irec),rdatahalf(1,1,k,ivar), &
                        nx,ny,header,maxlenh,rec_dathalf(k,ivar,irec),ierr_interp, &
                        iprint,bufix,bufjy,interp_type,label)
                   if (ierr_interp .ne. 0) nerr=nerr+1
                end do
             end do
          end if
          if (nvarfull .gt. 0) then
             do ivar=1,nvarfull
                do k=0,nlev_half
                   call interp_2d(rec_lon(irec),rec_lat(irec),rdatafull(1,1,k,ivar), &
                        nx,ny,header,maxlenh,rec_datfull(k,ivar,irec),ierr_interp, &
                        iprint,bufix,bufjy,interp_type,label)
                   if (ierr_interp .ne. 0) nerr=nerr+1
                end do
             end do
          end if
       end if
       if (rec_ind(irec) .eq. itime_data-1) then
          ! time level 2
          if (nvar2d .gt. 0) then
             do ivar=1,nvar2d
                call interp_2d(rec_lon(irec),rec_lat(irec),rdata2d(1,1,ivar),nx,ny, &
                     header,maxlenh,rec2_dat2d(ivar,irec),ierr_interp,iprint,bufix, &
                     bufjy,interp_type,label)
                if (ierr_interp .ne. 0) nerr=nerr+1
             end do
          end if
          if (nvarhalf .gt. 0) then
             do ivar=1,nvarhalf
                do k=1,nlev_half
                   call interp_2d(rec_lon(irec),rec_lat(irec),rdatahalf(1,1,k,ivar), &
                        nx,ny,header,maxlenh,rec2_dathalf(k,ivar,irec),ierr_interp, &
                        iprint,bufix,bufjy,interp_type,label)
                   if (ierr_interp .ne. 0) nerr=nerr+1
                end do
             end do
          end if
          if (nvarfull .gt. 0) then
             do ivar=1,nvarfull
                do k=0,nlev_half
                   call interp_2d(rec_lon(irec),rec_lat(irec),rdatafull(1,1,k,ivar), &
                        nx,ny,header,maxlenh,rec2_datfull(k,ivar,irec),ierr_interp, &
                        iprint,bufix,bufjy,interp_type,label)
                   if (ierr_interp .ne. 0) nerr=nerr+1
                end do
             end do
          end if
       end if
    end do

    if (nerr .ne. 0) then
       write (*,*) 'interp_data: internal logic error, nerr=',nerr
       stop 'interp_data: internal logic error'
    end if

  end subroutine interp_data

end PROGRAM interp_receptors
