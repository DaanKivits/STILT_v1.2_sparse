!-------------------------------------------------------------------------------
! Lists header contents for one netCDF file, and then print the data variable
! at position 1,1,1 at each time period.
!-------------------------------------------------------------------------------
! Last Revised: 06 Mar 2001
!-------------------------------------------------------------------------------

PROGRAM cdfhead

  IMPLICIT NONE

  REAL               :: level(20)
  INTEGER(8)         :: ctime
  INTEGER            :: n, k, kv, nz, ndat
  INTEGER            :: iy, im, id, ih, ic, julh
  INTEGER            :: fcopen , klen
  INTEGER(2)         :: kdata(144,73,17)
  REAL               :: rdata(144,73,17)
  CHARACTER(80)      :: fname

  CHARACTER(8)   :: varid
  INTEGER        :: handle      ! unit number
  INTEGER        :: kdate  (3)  ! starting year month day
  INTEGER        :: kdim   (4)  ! number in x,y,z,t
  INTEGER        :: vdata       ! byte position vertical levels
  INTEGER        :: nbyte       ! bytes per record variable    
  INTEGER        :: begin       ! byte position variable data   
  REAL           :: offset      ! packing offset value    
  REAL           :: scale       ! packing scaling value   

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
! open file for any variable

  varid='dummy&'
  WRITE(*,*)'Enter netCDF file name: '
  READ(5,'(A)')fname
  HANDLE=FCOPEN(fname,'r')

!---------------------------------------------------------
! scan each file for required starting information

  write(*,*)
  CALL ncscan (handle,varid,kdate,kdim,vdata,nbyte,begin,offset,scale,.true.)
  write(*,*)
  write(*,*)'Date: ',kdate
  write(*,*)'Dims: ',kdim

!---------------------------------------------------------
! construct time loop

  iy=kdate(1)
  im=kdate(2)
  id=kdate(3)
  ih=0
  ic=6
  
! initial position for all input data files
! back up one double integer for time record variaable
  CALL fcptps(handle,(begin-8),*900)

! last dimension is the number of time periods
  DO n=1,kdim(4)

!    load time variable
     CALL fcread(handle,ctime,8,1,*900)

!    load data (assume short integer)
     ndat=kdim(1)*kdim(2)*kdim(3)
     CALL fcread(handle,kdata,2,ndat,*900)

!    unpack the data from integer to real
     rdata =(kdata*scale)+offset

!    shift pointer
     begin=begin+nbyte

!    diagnostic output
     write(*,*)iy,im,id,ih,julh,' ... ',rdata(1,1,1) 
!    read(*,*)
     CALL NCJULH (iy,im,id,ih,ic,julh)

  END DO

!---------------------------------------------------------
! close all open files

  CALL FCCLOS(handle,*900)

  900 CONTINUE

END PROGRAM cdfhead 
