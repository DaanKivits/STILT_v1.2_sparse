! Use emacs f90-mode:   -*- mode:f90  -*-
module gettime_module

! $Id: gettime.f,v 1.2 2007/11/20 22:37:49 trn Exp $

contains
SUBROUTINE gettime (diag,ncid,itime,iyyyy,imo,ida,ihr,imn,iss,iyy,ntimes,ierr)

  use getdim_module
  use setvar_module

  IMPLICIT NONE

  include 'netcdf.inc'

  LOGICAL,       INTENT(IN)  :: diag
  INTEGER,       INTENT(IN)  :: ncid
  INTEGER,       INTENT(INOUT)  :: itime
  INTEGER,       INTENT(OUT) :: iyyyy,imo,ida,ihr,imn,iss,iyy,ntimes,ierr

  integer, parameter :: lchar=80
  integer :: DateStrLen
  CHARACTER(len=lchar) :: label
  CHARACTER(len=lchar) :: timestr
  CHARACTER(len=1),allocatable :: timevar(:,:)


  INTEGER       :: varid
  INTEGER       :: n3d  
  INTEGER       :: ndim
  INTEGER       :: dimlen(nf_max_var_dims)
  INTEGER       :: i,n,status
  INTEGER       :: dimids(nf_max_var_dims) 

  dimlen = 0
  varid  = 0
  n3d    = 0
  ndim   = 0
  ierr   = 0

  label = 'DateStrLen'
  call getdim(diag,ncid,label,DateStrLen,ierr)
  if (ierr .ne. 0) return
  if (DateStrLen .gt. lchar) then
     write (*,*) 'DateStrLen too long: ',DateStrLen
     ierr = 2
     return
  end if
  
  label = 'Times'
  call setvar(diag,ncid,label,varid,n3d,ndim,dimlen,ierr)
  ntimes = dimlen(ndim)
  if (ierr .ne. 0) return
  if (ndim .ne. 2 .or. dimlen(1) .ne. DateStrLen) &
       &  write (*,*) 'Unexpected dimensions for Times: ndim= ', &
       & ndim,' dimlen =',dimlen(1:ndim)
  if (itime .gt. ntimes) then
     write (*,*) 'Specified time period itime= ',itime, &
          & ' exceeds number of times in file, reset to ntimes=',ntimes
     itime=ntimes
  endif
  if (ndim .ne. 2 .or. dimlen(1) .ne. DateStrLen) then
     ierr = 3
     return
  end if
  
  allocate(timevar(dimlen(1),dimlen(2)),stat=status)
  if (status .ne. 0) then
     write (*,*) 'Error allocating timevar'
     ierr=4
     return
  end if

  status = NF_get_var_text(ncid, varid, timevar)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) 'Error getting timevar variable : ', NF_STRERROR(status)
     ierr=5
     return
  end IF
  timestr=' '
  do i=1,DateStrLen
     timestr(i:i) = timevar(i,itime)
  end do
  if (diag) write (*,*) 'Time string for period itime= ',itime,' : ',trim(timestr)

  ! decode string of the form: 2004-05-15_01:00:00
  read(timestr,'(i4,5(1x,i2))',iostat=status) iyyyy,imo,ida,ihr,imn,iss
  if (status .ne. 0) then
     write (*,*) 'Error decoding timestring:',trim(timestr)
     ierr=6
  end if
  iyy = mod(iyyyy,100)

  return

end SUBROUTINE gettime
end module gettime_module
