SUBROUTINE crs2dot (varcrs, vardot)

!-------------------------------------------------------------------------------
! Name:     Cross to Dot 
! Purpose:  Interpolates in horizontal from cross to dot points.
! Notes:    Based on PSU/NCAR model routine.
! Revised:  20 Apr 1999  Original version.  (TLO)
!           29 Oct 1999  Converted to free-form f90.  (TLO)
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER                      :: ix, jx, ie, je, i, j
  REAL,            INTENT(IN)  :: varcrs ( : , : )
  REAL,            INTENT(OUT) :: vardot ( : , : )

!-------------------------------------------------------------------------------
! Extract domain dimensions.
!-------------------------------------------------------------------------------

  IF ( SIZE(varcrs,1) /= SIZE(vardot,1) ) STOP 'crs2dot 1'
  IF ( SIZE(varcrs,2) /= SIZE(vardot,2) ) STOP 'crs2dot 2'

  ix = SIZE(varcrs,1)
  jx = SIZE(varcrs,2)

!!! new variables for test
  ie = ix - 1
  je = jx - 1

!-------------------------------------------------------------------------------
! For interior of grid, interpolate cross point values to dot points using
! four-point interpolation.
!-------------------------------------------------------------------------------

  DO j = 2, je
    DO i = 2, ie
      vardot(i,j) = 0.25 * ( varcrs(i,j)   + varcrs(i-1,j)  &
                           + varcrs(i,j-1) + varcrs(i-1,j-1) )
    ENDDO
  ENDDO

!!!  vardot(2:ix-1,2:jx-1) = 0.25 * ( varcrs(2:ix-1,2:jx-1) +  &
!!!                                   varcrs(1:ix-2,2:jx-1) +  &
!!!                                   varcrs(2:ix-1,1:jx-2) +  &
!!!                                   varcrs(1:ix-2,1:jx-2) )

!-------------------------------------------------------------------------------
! For outermost rows and columns, interpolate cross point values to dot points
! using two-point interpolation.  In row and column 1, there are no cross points
! below or to the left of the dow row that is being interpolated.  In row JX
! and column IX, cross points are not defined.
!-------------------------------------------------------------------------------

  DO i = 2, ie
    vardot(i,1)  = 0.5 * ( varcrs(i,1)  + varcrs(i-1,1)  )
    vardot(i,jx) = 0.5 * ( varcrs(i,je) + varcrs(i-1,je) )
  ENDDO

!!!  vardot(2:ix-1,1:jx:jx-1) = 0.5 * ( varcrs(2:ix-1,1:jx-1:jx-2) +  &
!!!                                     varcrs(1:ix-2,1:jx-1:jx-2) )

  DO j = 2, je
    vardot(1, j) = 0.5 * ( varcrs(1, j) + varcrs(1, j-1) )
    vardot(ix,j) = 0.5 * ( varcrs(ie,j) + varcrs(ie,j-1) )
  ENDDO

!!!  vardot(1:ix:ix-1,2:jx-1) = 0.5 * ( varcrs(1:ix-1:ix-2,2:jx-1) +  &
!!!                                     varcrs(1:ix-1:ix-2,1:jx-2) )

!-------------------------------------------------------------------------------
! Define dot point corners by persisting cross point corners.
!-------------------------------------------------------------------------------

  vardot(1, 1)  = varcrs(1, 1)
  vardot(1, jx) = varcrs(1, je)
  vardot(ix,jx) = varcrs(ie,je)
  vardot(ix,1)  = varcrs(ie,1)

!!!  vardot(1:ix:ix-1,1:jx:jx-1) = vardot(1:ix-1:ix-2,1:jx-1:jx-2)

  RETURN

END SUBROUTINE crs2dot
