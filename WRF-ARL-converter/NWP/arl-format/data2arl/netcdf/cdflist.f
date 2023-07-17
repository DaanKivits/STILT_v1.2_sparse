!-------------------------------------------------------------------------------
! Lists contents of all netCDF files that are required to create an ARL
! packed meteorological data set.
!-------------------------------------------------------------------------------
! Last Revised: 06 Mar 2001
!-------------------------------------------------------------------------------

PROGRAM cdflist

  IMPLICIT NONE

  INTEGER, PARAMETER :: nvar = 11       ! number of variables (and files)

  REAL               :: level(20)
  INTEGER(8)         :: ctime
  INTEGER            :: n, k, kv, nz, ndat
  INTEGER            :: iy, im, id, ih, ic, julh
  INTEGER            :: year, month, fcopen , klen
  INTEGER(2)         :: kdata(144,73,17)
  REAL               :: rdata(144,73,17)
  CHARACTER(80)      :: fname

!---------------------------------------------------------
! file dependent variables

  CHARACTER(8)   :: varb   (nvar)       ! data variable
  INTEGER        :: handle (nvar)       ! unit number
  INTEGER        :: kdate  (3,nvar)     ! starting year month day
  INTEGER        :: kdim   (4,nvar)     ! number in x,y,z,t
  INTEGER        :: vdata  (nvar)       ! byte position vertical levels
  INTEGER        :: nbyte  (nvar)       ! bytes per record variable    
  INTEGER        :: begin  (nvar)       ! byte position variable data   
  REAL           :: offset (nvar)       ! packing offset value    
  REAL           :: scale  (nvar)       ! packing scaling value   

!---------------------------------------------------------
! 5 surface variables followed by 6 3D variables

  DATA varb/'hgt&','pres&','air&','uwnd&','vwnd&',                &
            'hgt&','air&','uwnd&','vwnd&','omega&','rhum&'/

!---------------------------------------------------------
  INTERFACE
  SUBROUTINE ncscan (handle,varid,kdate,kdim,vdata,               &
                    nbyte,begin,offset,scale,diag)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: handle
  CHARACTER(8), INTENT(IN)  :: varid
  INTEGER,      INTENT(OUT) :: kdate (3)
  INTEGER,      INTENT(OUT) :: kdim  (4)
  INTEGER,      INTENT(OUT) :: vdata 
  INTEGER,      INTENT(OUT) :: nbyte 
  INTEGER,      INTENT(OUT) :: begin 
  REAL,         INTENT(OUT) :: offset
  REAL,         INTENT(OUT) :: scale
  LOGICAL,      INTENT(IN)  :: diag
  END SUBROUTINE ncscan 

  SUBROUTINE ncjulh (iy,im,id,ih,ic,julh)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: iy, im, id, ih
  INTEGER, INTENT(IN)    :: ic
  INTEGER, INTENT(OUT)   :: julh
  END SUBROUTINE ncjulh 
  END INTERFACE

!---------------------------------------------------------
! select the year and month for processing

  year=1987
  month=1 
        
!---------------------------------------------------------
! construct and open files based upon the year selected

  DO K=1,nvar
     klen=INDEX(varb(k),'&')-1
     IF(K.EQ.1)THEN
        WRITE(fname,'(a,a7)')varb(k)(1:klen),'.sfc.nc'
     ELSEIF(K.EQ.2)THEN
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.sfc.',year,'.nc'
     ELSEIF(K.GE.3.AND.K.LE.5)THEN
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.sig995.',year,'.nc'
     ELSE     
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.',year,'.nc'
     END IF
     WRITE(*,*)'Opening: ', fname(1:40)
     HANDLE(k)=FCOPEN(fname,'r')
  END DO

!---------------------------------------------------------
! scan each file for required starting information

  WRITE(*,*)
  DO k=1,nvar 
     CALL ncscan (handle(k),varb(k),kdate(:,k),kdim(:,k),vdata(k),             &
                 nbyte(k),begin(k),offset(k),scale(k),.false.)
     klen=INDEX(varb(k),'&')-1
     write(*,*)'Set: ',varb(k)(1:klen),kdate(:,k),kdim(:,k)
  END DO

!---------------------------------------------------------
! construct level information

  WRITE(*,*)
  WRITE(*,*)'Vertical data -'
  nz=0
  DO k=1,nvar
     IF(kdim(3,k).GT.nz)THEN
        nz=kdim(3,k)
        kv=k
     END IF
  END DO
  CALL fcptps(handle(kv),vdata(kv),*900)
  CALL fcread(handle(kv),level,4,nz,*900)
  k=nz
  DO WHILE (k.GT.0)
     write(*,*)k,' - ',level(k)
     k=k-1
  END DO

!---------------------------------------------------------
! construct time loop

  iy=kdate(1,kv)
  im=kdate(2,kv)
  id=kdate(3,kv)
  ih=0
  ic=6
  
! initial position for all input data files
  DO n=1,nvar 
!    back up one double integer for time record variaable
     CALL fcptps(handle(n),(begin(n)-8),*900)
  END DO

  WRITE(*,*)
  DO n=1,kdim(4,kv)
     WRITE(*,*)iy,im,id,ih

     vloop : DO k=1,nvar 
        IF(N.GT.1.AND.K.EQ.1) CYCLE vloop

!       load time variable
        CALL fcread(handle(k),ctime,8,1,*900)

!       load data (assume short integer)
        ndat=kdim(1,k)*kdim(2,k)*kdim(3,k)
        CALL fcread(handle(k),kdata,2,ndat,*900)

!       unpack the data from integer to real
        rdata =(kdata*scale(k))+offset(k)

!       shift pointer
        begin(k)=begin(k)+nbyte(k)

!       diagnostic output
        klen=INDEX(varb(k),'&')-1
        write(*,*)'    ',varb(k)(1:klen),ndat,nbyte(k), rdata(1,1,1) 
     END DO vloop

     CALL NCJULH (iy,im,id,ih,ic,julh)
  END DO

!---------------------------------------------------------
! close all open files

  DO K=1,nvar 
     CALL FCCLOS(handle(k),*900)
  END DO

  900 CONTINUE

END PROGRAM cdflist 
