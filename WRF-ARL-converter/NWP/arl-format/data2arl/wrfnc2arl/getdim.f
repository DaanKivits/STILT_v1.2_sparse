! Use emacs f90-mode:   -*- mode:f90  -*-
module getdim_module

! $Id: getdim.f,v 1.2 2007/06/19 15:44:16 trn Exp $

contains
SUBROUTINE getdim (diag,ncid,label,outlen,ierr)

  implicit none 

  include 'netcdf.inc'

  LOGICAL,       INTENT(IN)  :: diag
  INTEGER,       INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: label 
  INTEGER,       INTENT(OUT) :: outlen
  INTEGER,       INTENT(OUT) :: ierr

  integer :: status, dimid

  status = nf_inq_dimid(ncid,label,dimid)
  if (status .ne. nf_noerr) then
     write (*,*) 'Error getting dimid for',trim(label),' : ',nf_strerror(status)
     ierr = 1
  endif
  status = nf_inq_dimlen(ncid,dimid,outlen)
  if (status .ne. nf_noerr) then
     write (*,*) 'Error getting dimlen for',trim(label),' : ',nf_strerror(status)
     ierr = 2
  endif

  return
end SUBROUTINE getdim
end module getdim_module
