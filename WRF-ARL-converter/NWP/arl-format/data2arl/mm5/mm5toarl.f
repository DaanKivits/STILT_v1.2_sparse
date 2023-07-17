PROGRAM mm5toarl

!-------------------------------------------------------------------------------
! Name:     MM5 to ARL
! Purpose:  Get output from Penn State/NCAR Mesoscale Model (MM5) Version 3
!           and convert to ARL one-byte packed format.  
! Notes:    This program is based upon Tanya Otte's mm5tomlm.
! Revised:  04 May 2000  Original version.  (TLO)
!           14 Aug 2000  Converted to ARL   (RRD)
!           09 Mar 2004  Added TKE field    (RRD)
!-------------------------------------------------------------------------------

  USE file

  IMPLICIT NONE

  INTEGER,       PARAMETER     :: maxlevels  = 25    ! ARL format data output
  INTEGER,       PARAMETER     :: maxvariab  = 20    ! number different variab 

  CHARACTER*24                 :: currentdate
  REAL,          ALLOCATABLE   :: terrain    ( : , : )
  REAL,          ALLOCATABLE   :: halfsig    ( : )
  REAL,          ALLOCATABLE   :: arlsig     ( : )
  REAL,          ALLOCATABLE   :: hght       ( : , : , : )
  REAL,          ALLOCATABLE   :: hold       ( : , : )
  LOGICAL                      :: iend       = .FALSE.
  LOGICAL                      :: ifirst     = .TRUE.
  INTEGER                      :: ix
  INTEGER                      :: jx
  INTEGER                      :: kx
  INTEGER                      :: ka
  INTEGER                      :: n
  REAL                         :: p0s
  REAL,          ALLOCATABLE   :: pp         ( : , : , : )
  REAL,          ALLOCATABLE   :: precip     ( : , : )
  REAL,          ALLOCATABLE   :: pres       ( : , : , : )
  REAL,          ALLOCATABLE   :: psa        ( : , : )
  REAL,          ALLOCATABLE   :: psfc       ( : , : )
  REAL                         :: ptop
  REAL,          ALLOCATABLE   :: tke        ( : , : , : )
  REAL,          ALLOCATABLE   :: qva        ( : , : , : )
  REAL,          ALLOCATABLE   :: pblhgt     ( : , : )
  REAL,          ALLOCATABLE   :: shflux     ( : , : )
  REAL,          ALLOCATABLE   :: swdown     ( : , : )
  REAL,          ALLOCATABLE   :: ta         ( : , : , : )
  REAL                         :: time
  REAL                         :: tiso
  REAL                         :: tlp
  REAL                         :: ts0
  REAL,          ALLOCATABLE   :: tsfc       ( : , : )
  REAL,          ALLOCATABLE   :: ua         ( : , : , : )
  REAL,          ALLOCATABLE   :: ust        ( : , : )
  REAL,          ALLOCATABLE   :: va         ( : , : , : )
  REAL,          ALLOCATABLE   :: wa         ( : , : , : )

  INTEGER                      :: narg,iargc
  INTEGER                      :: nicg  
  INTEGER                      :: njcg 
  INTEGER                      :: proj
  INTEGER                      :: gnum
  REAL                         :: dist   
  REAL                         :: clat 
  REAL                         :: clon
  REAL                         :: tru1 
  REAL                         :: tru2 
  REAL                         :: pole 
  REAL                         :: cgi1 
  REAL                         :: cgj1 
  REAL                         :: cgi2 
  REAL                         :: cgj2 

  INTEGER                      :: n0,n1
  CHARACTER*4                  :: vchar0     (maxvariab)
  CHARACTER*4                  :: vchar1     (maxvariab)
  CHARACTER*80                 :: fname

!-------------------------------------------------------------------------------
! Configure subroutine interface argumment lists
!-------------------------------------------------------------------------------

  INTERFACE

    SUBROUTINE setupv3 (fname, ix, jx, kx, ptop, p0s, ts0, tlp, tiso,    &
                        nicg, njcg, proj, gnum, dist, clat, clon,        &
                        tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2) 
      IMPLICIT NONE
      CHARACTER*80,  INTENT(IN)    :: fname
      INTEGER,       INTENT(OUT)   :: ix
      INTEGER,       INTENT(OUT)   :: jx
      INTEGER,       INTENT(OUT)   :: kx
      REAL,          INTENT(OUT)   :: p0s
      REAL,          INTENT(OUT)   :: ptop
      REAL,          INTENT(OUT)   :: tiso
      REAL,          INTENT(OUT)   :: tlp
      REAL,          INTENT(OUT)   :: ts0
      INTEGER,       INTENT(OUT)   :: nicg  
      INTEGER,       INTENT(OUT)   :: njcg 
      INTEGER,       INTENT(OUT)   :: proj
      INTEGER,       INTENT(OUT)   :: gnum
      REAL,          INTENT(OUT)   :: dist   
      REAL,          INTENT(OUT)   :: clat 
      REAL,          INTENT(OUT)   :: clon
      REAL,          INTENT(OUT)   :: tru1 
      REAL,          INTENT(OUT)   :: tru2 
      REAL,          INTENT(OUT)   :: pole 
      REAL,          INTENT(OUT)   :: cgi1 
      REAL,          INTENT(OUT)   :: cgj1 
      REAL,          INTENT(OUT)   :: cgi2 
      REAL,          INTENT(OUT)   :: cgj2 
    END SUBROUTINE setupv3

    !-------------------------------------------------------------------------

    SUBROUTINE setgrid (ix, jx, kx, nicg, njcg, proj, gnum, dist, clat, clon, &
                        tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2,             &
                        ka, n0, n1, vchar0, vchar1, arlsig)
      IMPLICIT NONE
      INTEGER,       INTENT(IN)    :: ix
      INTEGER,       INTENT(IN)    :: jx
      INTEGER,       INTENT(IN)    :: kx
      INTEGER,       INTENT(IN)    :: nicg  
      INTEGER,       INTENT(IN)    :: njcg 
      INTEGER,       INTENT(IN)    :: proj
      INTEGER,       INTENT(IN)    :: gnum
      REAL,          INTENT(IN)    :: dist   
      REAL,          INTENT(IN)    :: clat 
      REAL,          INTENT(IN)    :: clon
      REAL,          INTENT(IN)    :: tru1 
      REAL,          INTENT(IN)    :: tru2 
      REAL,          INTENT(IN)    :: pole 
      REAL,          INTENT(IN)    :: cgi1 
      REAL,          INTENT(IN)    :: cgj1 
      REAL,          INTENT(IN)    :: cgi2 
      REAL,          INTENT(IN)    :: cgj2 
      INTEGER,       INTENT(IN)    :: ka,n0,n1
      CHARACTER*4,   INTENT(IN)    :: vchar0  ( : )
      CHARACTER*4,   INTENT(IN)    :: vchar1  ( : )
      REAL,          INTENT(IN)    :: arlsig  ( : )
    END SUBROUTINE setgrid

    !-------------------------------------------------------------------------

    SUBROUTINE rdmm5v3 (ix, jx, kx, ua, va, wa, ta, tke, qva, pp, psa,      &
                        shflux, ust, terrain, pblhgt, swdown, precip, hold, &
                        halfsig, currentdate, time, iend)
      IMPLICIT NONE
      CHARACTER*24,  INTENT(OUT)   :: currentdate
      REAL,          INTENT(OUT)   :: time
      REAL,          INTENT(OUT)   :: terrain    ( : , : )
      REAL,          INTENT(OUT)   :: halfsig    ( : )
      REAL,          INTENT(INOUT) :: hold       ( : , : )
      LOGICAL,       INTENT(INOUT) :: iend
      INTEGER,       INTENT(IN)    :: ix
      INTEGER,       INTENT(IN)    :: jx
      INTEGER,       INTENT(IN)    :: kx
      REAL,          INTENT(OUT)   :: pp         ( : , : , : )
      REAL,          INTENT(OUT)   :: precip     ( : , : )
      REAL,          INTENT(OUT)   :: psa        ( : , : )
      REAL,          INTENT(OUT)   :: tke        ( : , : , : )
      REAL,          INTENT(OUT)   :: qva        ( : , : , : )
      REAL,          INTENT(OUT)   :: pblhgt     ( : , : )
      REAL,          INTENT(OUT)   :: shflux     ( : , : )
      REAL,          INTENT(OUT)   :: swdown     ( : , : )
      REAL,          INTENT(OUT)   :: ta         ( : , : , : )
      REAL,          INTENT(OUT)   :: ua         ( : , : , : )
      REAL,          INTENT(OUT)   :: ust        ( : , : )
      REAL,          INTENT(OUT)   :: va         ( : , : , : )
      REAL,          INTENT(OUT)   :: wa         ( : , : , : )
    END SUBROUTINE rdmm5v3

    !-------------------------------------------------------------------------

    SUBROUTINE convert  (halfsig, hght, kx, pp, psa, pres, &
                         psfc, ptop, tlp, ts0, ta, wa, qva )
      IMPLICIT NONE
      REAL,          INTENT(IN)    :: halfsig    ( : )
      REAL,          INTENT(OUT)   :: hght       ( : , : , : )
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
    END SUBROUTINE convert 

    !-------------------------------------------------------------------------

    SUBROUTINE setvarb (ptop, n0, n1, vchar0, vchar1, halfsig, arlsig)
      IMPLICIT NONE
      REAL,         INTENT(IN)     :: ptop
      INTEGER,      INTENT(OUT)    :: n0,n1
      CHARACTER*4,  INTENT(OUT)    :: vchar0    ( : )
      CHARACTER*4,  INTENT(OUT)    :: vchar1    ( : )
      REAL,         INTENT(OUT)    :: halfsig   ( : )
      REAL,         INTENT(OUT)    :: arlsig    ( : )
    END SUBROUTINE setvarb

    !-------------------------------------------------------------------------

    SUBROUTINE wrtarl (ix, jx, kx, ka, n0, n1, vchar0, vchar1, psfc, &
               pblhgt, shflux, ust, swdown, precip, terrain, ua, va, &
               wa, ta, pres, tke, qva, currentdate, time, gnum )
      IMPLICIT NONE
      INTEGER,       INTENT(IN)    :: ix
      INTEGER,       INTENT(IN)    :: jx
      INTEGER,       INTENT(IN)    :: kx
      INTEGER,       INTENT(IN)    :: ka, n0, n1
      CHARACTER*4,   INTENT(IN)    :: vchar0     ( : )
      CHARACTER*4,   INTENT(IN)    :: vchar1     ( : )
      REAL,          INTENT(IN)    :: psfc       ( : , : )
      REAL,          INTENT(IN)    :: pblhgt     ( : , : )
      REAL,          INTENT(IN)    :: shflux     ( : , : )
      REAL,          INTENT(IN)    :: ust        ( : , : )
      REAL,          INTENT(IN)    :: swdown     ( : , : )
      REAL,          INTENT(IN)    :: precip     ( : , : )
      REAL,          INTENT(IN)    :: terrain    ( : , : )
      REAL,          INTENT(IN)    :: ua         ( : , : , : )
      REAL,          INTENT(IN)    :: va         ( : , : , : )
      REAL,          INTENT(IN)    :: wa         ( : , : , : )
      REAL,          INTENT(IN)    :: ta         ( : , : , : )
      REAL,          INTENT(IN)    :: pres       ( : , : , : )
      REAL,          INTENT(IN)    :: tke        ( : , : , : )
      REAL,          INTENT(IN)    :: qva        ( : , : , : )
      CHARACTER*24,  INTENT(IN)    :: currentdate
      REAL,          INTENT(IN)    :: time
      INTEGER,       INTENT(IN)    :: gnum
    END SUBROUTINE wrtarl
 END INTERFACE

!-------------------------------------------------------------------------------
! Check command line for input file name
!-------------------------------------------------------------------------------

  NARG=IARGC()
  IF(NARG.EQ.0)THEN
     fname='fort.10'
  ELSE
     narg=1
     CALL GETARG(narg,fname)
     IF(FNAME(1:2).EQ.'-h')THEN
        WRITE(*,*)'Usage: mm5toarl [mmout_domain{x} (fort.10)]'
        STOP
     END IF
  END IF

!-------------------------------------------------------------------------------
! Get domain attributes for MM5 domain.
!-------------------------------------------------------------------------------

  CALL setupv3 (fname, ix, jx, kx, ptop, p0s, ts0, tlp, tiso,    &
                nicg, njcg, proj, gnum, dist, clat, clon,        &
                tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2) 

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------

  ALLOCATE ( terrain (ix, jx)     )
  ALLOCATE ( halfsig (        kx) )
  ALLOCATE ( arlsig  (        kx) )
  ALLOCATE ( hold    (ix, jx)     )
  ALLOCATE ( hght    (ix, jx, kx) )
  ALLOCATE ( pp      (ix, jx, kx) )
  ALLOCATE ( precip  (ix, jx)     )
  ALLOCATE ( pres    (ix, jx, kx) )
  ALLOCATE ( psa     (ix, jx)     )
  ALLOCATE ( psfc    (ix, jx)     )
  ALLOCATE ( tke     (ix, jx, kx) )
  ALLOCATE ( qva     (ix, jx, kx) )
  ALLOCATE ( pblhgt  (ix, jx)     )
  ALLOCATE ( shflux  (ix, jx)     )
  ALLOCATE ( swdown  (ix, jx)     )
  ALLOCATE ( ta      (ix, jx, kx) )
  ALLOCATE ( tsfc    (ix, jx)     )
  ALLOCATE ( ua      (ix, jx, kx) )
  ALLOCATE ( ust     (ix, jx)     )
  ALLOCATE ( va      (ix, jx, kx) )
  ALLOCATE ( wa      (ix, jx, kx+1))

!-------------------------------------------------------------------------------
! Set up ARL format variables and grid definitions
!-------------------------------------------------------------------------------

  CALL tminit
  CALL setvarb (ptop, n0, n1, vchar0, vchar1, halfsig, arlsig)

  ka = MIN ( maxlevels, kx )

  CALL setgrid (ix, jx, kx, nicg, njcg, proj, gnum, dist, clat, clon, &
                tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2,             &
                ka, n0, n1, vchar0, vchar1, arlsig)

  hold  = 0.0

!-------------------------------------------------------------------------------
! Loop over times in MM5 file.
!-------------------------------------------------------------------------------

  dataloop: DO

    CALL rdmm5v3 (ix, jx, kx, ua, va, wa, ta, tke, qva, pp, psa,   &
                  shflux, ust, terrain, pblhgt, swdown,            &
                  precip, hold, halfsig, currentdate, time, iend)

    IF ( iend ) EXIT dataloop

    CALL convert  (halfsig, hght, kx, pp, psa, pres, psfc, ptop,   &
                   tlp, ts0, ta, wa, qva )

    CALL wrtarl (ix, jx, kx, ka, n0, n1, vchar0, vchar1, psfc,     &
         pblhgt, shflux, ust, swdown, precip, terrain, ua, va,     &
         wa, ta, pres, tke, qva, currentdate, time, gnum )

  ENDDO dataloop

!-------------------------------------------------------------------------------
! Deallocate variables.
!-------------------------------------------------------------------------------

  DEALLOCATE ( terrain )
  DEALLOCATE ( halfsig )
  DEALLOCATE ( arlsig )
  DEALLOCATE ( hold    )
  DEALLOCATE ( hght    )
  DEALLOCATE ( pp      )
  DEALLOCATE ( precip  )
  DEALLOCATE ( psfc    )
  DEALLOCATE ( pres    )
  DEALLOCATE ( psa     )
  DEALLOCATE ( tke     )
  DEALLOCATE ( qva     )
  DEALLOCATE ( pblhgt  )
  DEALLOCATE ( shflux  )
  DEALLOCATE ( swdown  )
  DEALLOCATE ( ta      )
  DEALLOCATE ( tsfc    )
  DEALLOCATE ( ua      )
  DEALLOCATE ( ust     )
  DEALLOCATE ( va      )
  DEALLOCATE ( wa      )
 
  CLOSE (iutarl)
  CLOSE (iutmm )

END PROGRAM mm5toarl
