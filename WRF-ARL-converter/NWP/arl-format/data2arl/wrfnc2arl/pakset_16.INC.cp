INTERFACE 

! interface block for pakset_16 subroutine
! needed because we are using an optional argument

SUBROUTINE PAKSET_16(LUNIT,FNAME,KREC1,NXP,NYP,NZP,use_16,levfmt)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,       INTENT(IN)    :: lunit     ! output unit number
  CHARACTER(*),  INTENT(INOUT) :: fname     ! file name of METDATA.CFG
  INTEGER,       INTENT(IN)    :: krec1     ! position of index record at time-1
  INTEGER,       INTENT(OUT)   :: nxp, nyp  ! horizontal grid dimensions
  INTEGER,       INTENT(OUT)   :: nzp       ! vertical grid dimension (incl sfc)
  logical, intent (in)         :: use_16    ! sixteen bit file (double the precision)
  character (len=*), intent(in), optional :: levfmt

!-------------------------------------------------------------------------------

END SUBROUTINE pakset_16

end INTERFACE
