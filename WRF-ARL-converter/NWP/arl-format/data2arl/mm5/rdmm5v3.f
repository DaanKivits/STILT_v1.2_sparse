SUBROUTINE rdmm5v3 (ix, jx, kx, ua, va, wa, ta, tke, qva, pp, psa,  &
                    shflux, ust, terrain, pblhgt, swdown,           &
                    precip, hold, halfsig, currentdate, time, iend)

!-------------------------------------------------------------------------------
! Name:    Read MM5v3 Output
! Purpose: Reads MM5v3 output.
! Revised: 04 May 2000  Original version.  (TLO)
!          09 Aug 2000  conversion for ARL (RRD)
!          05 Sep 2000  Added cross to dot point interpolation. (RRD)
!          09 Mar 2004  Added TKE field (RRD)
!-------------------------------------------------------------------------------

  USE file

  IMPLICIT NONE

  INTEGER,       ALLOCATABLE   :: bhi       ( : , : )
  REAL,          ALLOCATABLE   :: bhr       ( : , : )
  CHARACTER*24,  INTENT(OUT)   :: currentdate
  REAL,          ALLOCATABLE   :: data      ( : , : , : , : )
  CHARACTER*46                 :: description
  INTEGER                      :: end_index ( 4 )
  LOGICAL                      :: gotterrain
  LOGICAL                      :: gotpp
  LOGICAL                      :: gotpsa
  LOGICAL                      :: gottke
  LOGICAL                      :: gotqva
  LOGICAL                      :: gotrainc
  LOGICAL                      :: gotrainn
  LOGICAL                      :: gotpblhgt
  LOGICAL                      :: gotshflux
  LOGICAL                      :: gotsig
  LOGICAL                      :: gotswdown
  LOGICAL                      :: gotta
  LOGICAL                      :: gotua
  LOGICAL                      :: gotust
  LOGICAL                      :: gotva
  LOGICAL                      :: gotwa  
  REAL,          INTENT(OUT)   :: halfsig   ( : )
  CHARACTER*24                 :: hdate
  REAL,          INTENT(INOUT) :: hold      ( : , : )
  INTEGER                      :: i
  INTEGER                      :: idts_end
  INTEGER                      :: idts_start
  LOGICAL,       INTENT(INOUT) :: iend
  LOGICAL,       SAVE          :: ifirst    = .TRUE.
  INTEGER                      :: iflag
  INTEGER                      :: istat
  INTEGER,       INTENT(IN)    :: ix
  INTEGER                      :: j
  INTEGER,       INTENT(IN)    :: jx
  INTEGER                      :: k
  INTEGER,       INTENT(IN)    :: kx
  CHARACTER*9                  :: name
  INTEGER                      :: ndim
  LOGICAL                      :: newtime
  INTEGER,       PARAMETER     :: numprogs  = 20
  INTEGER,       PARAMETER     :: numvalsi  = 50
  INTEGER,       PARAMETER     :: numvalsr  = 20
  CHARACTER*4                  :: ordering
  REAL,          INTENT(OUT)   :: pblhgt    ( : , : )
  REAL,          INTENT(OUT)   :: pp        ( : , : , : )
  REAL,          INTENT(OUT)   :: precip    ( : , : )
  REAL,          INTENT(OUT)   :: psa       ( : , : )
  REAL,          INTENT(OUT)   :: tke       ( : , : , : )
  REAL,          INTENT(OUT)   :: qva       ( : , : , : )
  REAL,          ALLOCATABLE   :: raincon   ( : , : )
  REAL,          ALLOCATABLE   :: rainnon   ( : , : )
  INTEGER,       SAVE          :: season
  REAL,          INTENT(OUT)   :: shflux    ( : , : )
  CHARACTER*4                  :: staggering
  INTEGER                      :: start_index(4)
  CHARACTER*24,  SAVE          :: startdate
  REAL,          INTENT(OUT)   :: swdown    ( : , : )
  REAL,          INTENT(OUT)   :: ta        ( : , : , : )
  REAL,          INTENT(OUT)   :: terrain   ( : , : )
  REAL,          INTENT(OUT)   :: time
  REAL,          INTENT(OUT)   :: ua        ( : , : , : )
  CHARACTER*25                 :: units
  REAL,          INTENT(OUT)   :: ust       ( : , : )
  REAL,          INTENT(OUT)   :: va        ( : , : , : )
  REAL,          INTENT(OUT)   :: wa        ( : , : , : )
  CHARACTER*9                  :: var

  INTERFACE

    SUBROUTINE crs2dot (varcrs, vardot)
      IMPLICIT NONE
      REAL,          INTENT(IN)    :: varcrs     ( : , : )
      REAL,          INTENT(OUT)   :: vardot     ( : , : )
    END SUBROUTINE crs2dot

  END INTERFACE

!-------------------------------------------------------------------------------
! Allocate necessary arrays.
!-------------------------------------------------------------------------------

  ALLOCATE ( bhi (numvalsi, numprogs) )
  ALLOCATE ( bhr (numvalsr, numprogs) )

  ALLOCATE ( raincon (ix, jx) )
  ALLOCATE ( rainnon (ix, jx) )

!-------------------------------------------------------------------------------
! Initialize data capture flags.  Need to do this here (rather than on
! declaration line) since Sun seems to default all variables to "SAVE".
! That is, on Sun, the variables will not be re-initialized from the
! declaration line on subsequent calls; they hold the value from the last call.
!-------------------------------------------------------------------------------

  gotterrain = .FALSE.
  gotpp      = .FALSE.
  gotpsa     = .FALSE.
  gottke     = .FALSE.
  gotqva     = .FALSE.
  gotrainc   = .FALSE.
  gotrainn   = .FALSE.
  gotpblhgt  = .FALSE.
  gotshflux  = .FALSE.
  gotsig     = .FALSE.
  gotswdown  = .FALSE.
  gotta      = .FALSE.
  gotua      = .FALSE.
  gotust     = .FALSE.
  gotva      = .FALSE.
  gotwa      = .FALSE.
  newtime    = .TRUE.

!-------------------------------------------------------------------------------
! Read MM5v3 data for this domain.  Note that big header character information
! (BHIC and BHRC) is not read and not used here. Horizontal wind velocity values
! remain on dot points, while thermodynamic variables interpolated from cross
! points to dot points.     
!-------------------------------------------------------------------------------

  v3data: DO
    var = 'IFLAG    '
    READ (iutmm, IOSTAT=istat, ERR=8000, END=999) iflag

    IF ( iflag == 0 ) THEN

      var = 'BHI & BHR'
      READ (iutmm, IOSTAT=istat, ERR=8000, END=8100) bhi, bhr

    ELSE IF ( iflag == 1 ) THEN

      var = 'SM HEADER'
      READ (iutmm, IOSTAT=istat, ERR=8000, END=8200) ndim, start_index,  &
            end_index, time, staggering, ordering, currentdate, name,    &
            units, description

      IF ( ifirst ) THEN
        startdate = currentdate
        ifirst = .FALSE.
      ENDIF

      IF ( newtime ) THEN
        WRITE (*,'(/,a,2x,f15.5," Hours"/)') currentdate, time/60.0
        newtime = .FALSE.
      ENDIF

      IF ( ndim == 1 ) THEN
        ALLOCATE ( data(end_index(1), 1, 1, 1) )
      ELSE IF ( ndim == 2 ) THEN
        ALLOCATE ( data(end_index(1), end_index(2), 1, 1) )
      ELSE IF ( ndim == 3 ) THEN
        ALLOCATE ( data(end_index(1), end_index(2), end_index(3), 1) )
      ELSE IF ( ndim == 4 ) THEN
        ALLOCATE ( data(end_index(1), end_index(2), end_index(3), end_index(4)))
      ENDIF

      var = name
      READ (iutmm, IOSTAT=istat, ERR=8000, END=8300) data

! u velocity component (m/s)

      IF ( ( name == 'U        ' ) .AND. ( .NOT. gotua ) ) THEN
        IF ( ( SIZE(ua,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(ua,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(ua,3) == SIZE(data,3) ) ) THEN
          ua = data(:,:,:,1)
          gotua = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! v velocity component (m/s)

      IF ( ( name == 'V        ' ) .AND. ( .NOT. gotva ) ) THEN
        IF ( ( SIZE(va,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(va,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(va,3) == SIZE(data,3) ) ) THEN
          va = data(:,:,:,1)
          gotva = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! w velocity component (m/s)

      IF ( ( name == 'W        ' ) .AND. ( .NOT. gotwa ) ) THEN
        IF ( ( SIZE(wa,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(wa,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(wa,3) == SIZE(data,3) ) ) THEN
          DO k = 1 , kx+1
             CALL crs2dot (data(:,:,k,1),wa(:,:,k))
          END DO
          gotwa = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! air temperature (deg K)

      IF ( ( name == 'T        ' ) .AND. ( .NOT. gotta ) ) THEN
        IF ( ( SIZE(ta,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(ta,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(ta,3) == SIZE(data,3) ) ) THEN
          DO k = 1 , kx 
             CALL crs2dot (data(:,:,k,1),ta(:,:,k))
          END DO
          gotta = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! turbulent kinetic energy (J/kg)

      IF ( ( name == 'TKE      ' ) .AND. ( .NOT. gottke ) ) THEN
        IF ( ( SIZE(tke,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(tke,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(tke,3) == SIZE(data,3) ) ) THEN
          DO k = 1 , kx 
             CALL crs2dot (data(:,:,k,1),tke(:,:,k))
          END DO
          gottke = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! water vapor mixing ratio (kg/kg)

      IF ( ( name == 'Q        ' ) .AND. ( .NOT. gotqva ) ) THEN
        IF ( ( SIZE(qva,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(qva,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(qva,3) == SIZE(data,3) ) ) THEN
          DO k = 1 , kx 
             CALL crs2dot (data(:,:,k,1),qva(:,:,k))
          END DO
          gotqva = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! pressure perturbation (Pa)

      IF ( ( name == 'PP       ' ) .AND. ( .NOT. gotpp ) ) THEN
        IF ( ( SIZE(pp,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(pp,2) == SIZE(data,2) ) .AND.  &
             ( SIZE(pp,3) == SIZE(data,3) ) ) THEN
          DO k = 1 , kx 
             CALL crs2dot (data(:,:,k,1),pp(:,:,k))
          END DO
          gotpp = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! reference surface pressure minus ptop (Pa)

      IF ( ( name == 'PSTARCRS ' ) .AND. ( .NOT. gotpsa ) ) THEN
        IF ( ( SIZE(psa,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(psa,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),psa)
          gotpsa = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF
 
! accumulated convective precipitation (cm)

      IF ( ( name == 'RAIN CON ' ) .AND. ( .NOT. gotrainc ) ) THEN
        IF ( ( SIZE(raincon,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(raincon,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),raincon)
          raincon = 0.01*raincon                ! convert to meters
          gotrainc = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! accumulated non-convective precipitation (cm)

      IF ( ( name == 'RAIN NON ' ) .AND. ( .NOT. gotrainn ) ) THEN
        IF ( ( SIZE(rainnon,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(rainnon,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),rainnon)
          rainnon = 0.01*rainnon                ! convert to meters
          gotrainn = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! sensible heat flux (W/m^2)

      IF ( ( name == 'SHFLUX   ' ) .AND. ( .NOT. gotshflux ) ) THEN
        IF ( ( SIZE(shflux,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(shflux,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),shflux)
          gotshflux = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! friction velocity (m/s)

      IF ( ( name == 'UST      ' ) .AND. ( .NOT. gotust ) ) THEN
        IF ( ( SIZE(ust,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(ust,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),ust)
          gotust = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! downward shortwave radiation (W/m^2)

      IF ( ( name == 'SWDOWN   ' ) .AND. ( .NOT. gotswdown ) ) THEN
        IF ( ( SIZE(swdown,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(swdown,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),swdown)
          gotswdown = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! terrain height (meters)

      IF ( ( name == 'TERRAIN  ' ) .AND. ( .NOT. gotterrain ) ) THEN
        IF ( ( SIZE(terrain,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(terrain,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),terrain)
          gotterrain = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! pbl height (meters)

      IF ( ( name == 'PBL HGT  ' ) .AND. ( .NOT. gotpblhgt ) ) THEN
        IF ( ( SIZE(pblhgt,1) == SIZE(data,1) ) .AND.  &
             ( SIZE(pblhgt,2) == SIZE(data,2) ) ) THEN
          CALL crs2dot (data(:,:,1,1),pblhgt)
          gotpblhgt = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

! sigma profile (index = 1 at top)

      IF ( ( name == 'SIGMAH   ' ) .AND. ( .NOT. gotsig ) ) THEN
        IF ( SIZE(halfsig) == SIZE(data,1) ) THEN
          halfsig = data(:,1,1,1)
          gotsig = .TRUE.
        ELSE
          GOTO 8400
        ENDIF
      ENDIF

      DEALLOCATE ( data )

    ELSE IF ( iflag == 2 ) THEN

      newtime = .TRUE.
      EXIT v3data   ! only need one time period in this subroutine call

    ELSE

      PRINT*, '*** UNKNOWN MM5v3 DATA FLAG...IFLAG = ', iflag
      STOP

    ENDIF

    CYCLE v3data
 999 iend = .TRUE.
     EXIT v3data

  ENDDO v3data

!-------------------------------------------------------------------------------
! Make sure we collected the arrays we need.
! Missing TKE is not a failure ... may not be in original file
!-------------------------------------------------------------------------------

  IF ( iend ) RETURN

  IF ( .NOT. gotua ) THEN
    PRINT*, 'DID NOT GET U ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotva ) THEN
    PRINT*, 'DID NOT GET V ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotwa ) THEN
    PRINT*, 'DID NOT GET W ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotta ) THEN
    PRINT*, 'DID NOT GET T ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotqva ) THEN
    PRINT*, 'DID NOT GET Q ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotpp ) THEN
    PRINT*, 'DID NOT GET PP ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotpsa ) THEN
    PRINT*, 'DID NOT GET PSTARCRS ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotrainc ) THEN
    PRINT*, 'DID NOT GET RAIN CON ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotrainn ) THEN
    PRINT*, 'DID NOT GET RAIN NON ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotshflux ) THEN
    PRINT*, 'DID NOT GET SHFLUX ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotust ) THEN
    PRINT*, 'DID NOT GET UST ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotswdown ) THEN
    PRINT*, 'DID NOT GET SWDOWN ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotterrain ) THEN
    PRINT*, 'DID NOT GET TERRAIN ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotpblhgt ) THEN
    PRINT*, 'DID NOT GET PBL HGT ARRAY'
    STOP
  ENDIF

  IF ( .NOT. gotsig ) THEN
    PRINT*, 'DID NOT GET SIGMAH ARRAY'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Fill precipitation array (PRECIP) with amount over time interval.
!-------------------------------------------------------------------------------

  precip = raincon + rainnon - hold  ! "bucket" total in cm

  hold   = raincon + rainnon         ! accumulated total to use next interval

!-------------------------------------------------------------------------------
! Convert pressure variables from Pascals (Pa) to millibars (hPa).
!-------------------------------------------------------------------------------

  pp  = pp  / 100.0
  psa = psa / 100.0

!-------------------------------------------------------------------------------
! Deallocate arrays.
!-------------------------------------------------------------------------------

  DEALLOCATE ( bhi     )
  DEALLOCATE ( bhr     )
  DEALLOCATE ( raincon )
  DEALLOCATE ( rainnon )

  RETURN

!-------------------------------------------------------------------------------
! Error-handling section.
!-------------------------------------------------------------------------------

 8000 WRITE (6,9000) iutmm, TRIM(var), istat
      CALL abort
      RETURN

 8100 WRITE (6,9100) iutmm, iflag, istat
      CALL abort
      RETURN

 8200 WRITE (6,9200) iutmm, istat
      CALL abort
      RETURN

 8300 WRITE (6,9300) iutmm, TRIM(var), istat
      CALL abort
      RETURN

 8400 WRITE (6,9400) TRIM(var)
      CALL abort
      RETURN

 9000 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** SUBROUTINE: RDMM5V3',                   &
              /, 1x, '***   ERROR READING FILE, UNIT = ', i3,     &
              /, 1x, '***   VARIABLE = ', a,                      &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9100 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** SUBROUTINE: RDMM5V3',                   &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   IFLAG = ', i3,                        &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9200 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** SUBROUTINE: RDMM5V3',                   &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   VARIABLE = SMALL HEADER',             &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9300 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** SUBROUTINE: RDMM5V3',                   &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   VARIABLE = ', a,                      &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9400 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** SUBROUTINE: RDMM5V3',                   &
              /, 1x, '***   FOUND VARIABLE ', a,                  &
              /, 1x, '***   BUT ARRAY DIMENSIONS DO NOT MATCH',   &
              /, 1X, 70('*'))

END SUBROUTINE rdmm5v3
