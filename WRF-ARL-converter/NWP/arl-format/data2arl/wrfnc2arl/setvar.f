! Use emacs f90-mode:   -*- mode:f90  -*-
module setvar_module

! $Id: setvar.f,v 1.2 2007/06/19 15:44:16 trn Exp $

contains
SUBROUTINE setvar (diag,ncid,label,varid,n3d,ndim,dimlen,ierr)

  IMPLICIT NONE

  include 'netcdf.inc'

  LOGICAL,       INTENT(IN)  :: diag
  INTEGER,       INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: label 
  INTEGER,       INTENT(OUT) :: varid
  INTEGER,       INTENT(OUT) :: n3d  
  INTEGER,       INTENT(OUT) :: ndim
  INTEGER,       INTENT(OUT) :: dimlen(nf_max_var_dims)
  INTEGER,       INTENT(OUT) :: ierr

  INTEGER       :: n,status
  INTEGER       :: dimids(nf_max_var_dims) 

  dimlen = 0
  varid  = 0
  n3d    = 0
  ndim   = 0
  ierr   = 0

  status = NF_INQ_VARID(ncid, label, varid)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),'= ',NF_STRERROR(status)
     ierr=1
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),'= ',varid
  END IF

  status = NF_INQ_VARNDIMS(ncid, varid, ndim)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),'= ',NF_STRERROR(status)
     ierr=2
     return
  ELSE
     IF(diag)WRITE(*,*)' Number of dimensions = ',ndim
  END IF

  status = NF_INQ_VARDIMID(ncid, varid, dimids)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),'= ',NF_STRERROR(status)
     ierr=3
     return
  ELSE
     IF(diag)WRITE(*,*)'Dimension IDs = ',dimids(1:ndim)
  END IF

  DO n=1,ndim
     status = NF_INQ_DIMLEN(ncid, dimids(n), dimlen(n))
     IF(status.NE.nf_noerr) THEN
        WRITE(*,*) trim(label),'= ',NF_STRERROR(status)
        ierr=4
        return
     ELSE
        IF(diag)WRITE(*,*)'n = ',n,' Length = ',dimlen(n)
     END IF
  END DO

! additional information for 3D variables
  IF(ndim.GT.2)THEN
     status = NF_GET_ATT_INT(ncid, varid, '_n3D', n3d)
     IF(status.NE.nf_noerr) THEN
!###    WRITE(*,*) '_n3D: ',NF_STRERROR(status)
        n3d=0     
     ELSE
        IF(diag)WRITE(*,*)'_n3D = ',n3d
     END IF

  END IF

END SUBROUTINE setvar 
end module setvar_module
