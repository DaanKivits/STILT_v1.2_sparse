SUBROUTINE convert  (halfsig, hght, kx, pp, psa, pres, psfc, ptop, &
                     tlp, ts0, ta, wa, qva)

!-------------------------------------------------------------------------------
! Name:     Convert      
! Purpose:  Compute surface and pressure profile from perturbation pressure and
!           converts vertical velocity (m/s) to omega units (hPa/s).
! Notes:    Based upon equations in TLO's sfclayer subroutine 
! Revised:  14 Aug 2000  Original version.  (RRD)
!           05 Sep 2000  Vertical velocity offset by one index. (RRD)
!           03 Sep 2002  Pressure computation correction (KAG)
!-------------------------------------------------------------------------------

  USE consts

  IMPLICIT NONE

  REAL,          INTENT(IN)    :: halfsig    ( : )
  REAL,          INTENT(OUT)   :: hght       ( : , : , : )
  INTEGER                      :: k
  INTEGER,       INTENT(IN)    :: kx
  REAL,          INTENT(IN)    :: pp         ( : , : , : )
  REAL,          INTENT(IN)    :: psa        ( : , : )
  REAL,          INTENT(OUT)   :: pres       ( : , : , : )
  REAL,          INTENT(OUT)   :: psfc       ( : , : )
  REAL,          INTENT(IN)    :: ptop
  REAL,          INTENT(IN)    :: tlp
  REAL,          INTENT(IN)    :: ts0
  REAL,          INTENT(IN)    :: ta         ( : , : , : )
  REAL,          INTENT(INOUT) :: wa         ( : , : , : )
  REAL,          INTENT(IN)    :: qva        ( : , : , : )

!-------------------------------------------------------------------------------
! Pressure and Height Computation
!-------------------------------------------------------------------------------

  psfc(:,:) =  psa(:,:)  + ptop + pp(:,:,kx) 

  DO k=1,kx

     pres(:,:,k) = (halfsig(k) * psa(:,:)) + ptop + pp(:,:,k)

     hght(:,:,k) = - (r * tlp / (2.0 * g) * (ALOG(pres(:,:,k)/psfc(:,:))) ** 2 &
                   +  r * ts0 /        g  *  ALOG(pres(:,:,k)/psfc(:,:)))

  ! convert vertical velocity from m/s to hPa/s using omega = -W g rho

     wa(:,:,k) = -wa(:,:,k+1) * g * pres(:,:,k) /  &
                 (r * ta(:,:,k) * ( 1.0 + 0.6077 * qva(:,:,k) ) )

  END DO

  RETURN

END SUBROUTINE convert 
