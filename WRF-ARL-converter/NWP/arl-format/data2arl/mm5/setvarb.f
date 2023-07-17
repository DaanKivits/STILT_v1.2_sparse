SUBROUTINE setvarb (ptop, n0, n1, vchar0, vchar1, halfsig, arlsig)

!-------------------------------------------------------------------------------
! Name:     Setvarb       
! Purpose:  Sets variable arrays for ARL packing format based upon V3 contents.
! Notes:    Reads variables from the first time period then rewinds.         
! Revised:  14 Aug 2000 (RRD) - Original version
!           09 Mar 2004 (RRD) - Added TKE as output field 
!-------------------------------------------------------------------------------

  USE file

  IMPLICIT NONE

  INTEGER,       ALLOCATABLE   :: bhi       ( : , : )
  CHARACTER*80,  ALLOCATABLE   :: bhic      ( : , : )
  REAL,          ALLOCATABLE   :: bhr       ( : , : )
  CHARACTER*80,  ALLOCATABLE   :: bhrc      ( : , : )
  CHARACTER*24                 :: currentdate
  REAL,          ALLOCATABLE   :: data      ( : , : , : , : )
  CHARACTER*46                 :: description
  INTEGER                      :: end_index ( 4 )
  INTEGER                      :: iflag
  INTEGER                      :: istat
  CHARACTER*9                  :: name, var
  INTEGER                      :: ndim
  INTEGER,       PARAMETER     :: numprogs  = 20
  INTEGER,       PARAMETER     :: numvalsi  = 50
  INTEGER,       PARAMETER     :: numvalsr  = 20
  CHARACTER*4                  :: ordering
  REAL                         :: sample
  CHARACTER*4                  :: staggering
  INTEGER                      :: start_index ( 4 )
  REAL                         :: time
  CHARACTER*25                 :: units

  REAL,         INTENT(IN)     :: ptop
  INTEGER,      INTENT(OUT)    :: n0,n1
  CHARACTER*4,  INTENT(OUT)    :: vchar0    ( : )
  CHARACTER*4,  INTENT(OUT)    :: vchar1    ( : )
  REAL,         INTENT(OUT)    :: halfsig   ( : )
  REAL,         INTENT(OUT)    :: arlsig    ( : )

!-------------------------------------------------------------------------------
! Set variable availabity flags
!-------------------------------------------------------------------------------

  LOGICAL                      :: gotua
  LOGICAL                      :: gotva
  LOGICAL                      :: gotwa
  LOGICAL                      :: gotta
  LOGICAL                      :: gottke
  LOGICAL                      :: gotqva
  LOGICAL                      :: gotpp
  LOGICAL                      :: gotrain 
  LOGICAL                      :: gotshflux
  LOGICAL                      :: gotust
  LOGICAL                      :: gotswdown
  LOGICAL                      :: gotterrain
  LOGICAL                      :: gotpblhgt
  LOGICAL                      :: gotsig

  gotua      = .FALSE.
  gotva      = .FALSE.
  gotwa      = .FALSE.
  gotta      = .FALSE.
  gottke     = .FALSE.
  gotqva     = .FALSE.
  gotpp      = .FALSE.
  gotrain    = .FALSE.
  gotshflux  = .FALSE.
  gotust     = .FALSE.
  gotswdown  = .FALSE.
  gotterrain = .FALSE.
  gotpblhgt  = .FALSE.
  gotsig     = .FALSE.

!-------------------------------------------------------------------------------
! Allocate header arrays.
!-------------------------------------------------------------------------------

  ALLOCATE ( bhi  (numvalsi, numprogs) )
  ALLOCATE ( bhic (numvalsi, numprogs) )
  ALLOCATE ( bhr  (numvalsr, numprogs) )
  ALLOCATE ( bhrc (numvalsr, numprogs) )

!-------------------------------------------------------------------------------
! Loop through the first time period in the MM5v3 output file.
!-------------------------------------------------------------------------------

  v3data: DO

    var = 'IFLAG    '
    READ (iutmm, IOSTAT=istat, ERR=8000, END=8300) iflag

    IF ( iflag == 0 ) THEN

      var = 'BIG HEADR'
      READ (iutmm, IOSTAT=istat, ERR=8000, END=8100) bhi, bhr, bhic, bhrc

      n0 = 0  ! surface variable counter
      n1 = 0  ! upper level variable counter

    ELSE IF ( iflag == 1 ) THEN

      var = 'SM HEADER'
      READ (iutmm, IOSTAT=istat, ERR=8000, END=8200) ndim, start_index,  &
            end_index, time, staggering, ordering, currentdate, name,    &
            units, description

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

      IF ( ( name == 'U        ' ) .AND. ( .NOT. gotua ) ) THEN
          gotua = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'UWND'
      END IF

      IF ( ( name == 'V        ' ) .AND. ( .NOT. gotva ) ) THEN
          gotva = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'VWND'
      END IF

      IF ( ( name == 'W        ' ) .AND. ( .NOT. gotwa ) ) THEN
          gotva = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'WWND'
      END IF

      IF ( ( name == 'T        ' ) .AND. ( .NOT. gotta ) ) THEN
          gotta = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'TEMP'
      END IF

      IF ( ( name == 'Q        ' ) .AND. ( .NOT. gotqva ) ) THEN
          gotqva = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'SPHU'
      END IF

      IF ( ( name == 'TKE      ' ) .AND. ( .NOT. gottke ) ) THEN
          gottke = .TRUE.
          n1 = n1 + 1
          vchar1 (n1) = 'TKEN'
      END IF

      IF ( ( name == 'PP       ' ) .AND. ( .NOT. gotpp ) ) THEN
          gotpp = .TRUE.
          n1 = n1 + 1
          n0 = n0 + 1
          vchar0 (n0) = 'PRSS'
          vchar1 (n1) = 'PRES'
      END IF

      IF ( ( name(1:4) == 'RAIN' ) .AND. ( .NOT. gotrain ) ) THEN
          gotrain = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'TPP1'
      END IF

      IF ( ( name == 'SHFLUX   ' ) .AND. ( .NOT. gotshflux ) ) THEN
          gotshflux = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'SHTF'
      END IF

      IF ( ( name == 'UST      ' ) .AND. ( .NOT. gotust ) ) THEN
          gotust = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'USTR'
      END IF

      IF ( ( name == 'SWDOWN   ' ) .AND. ( .NOT. gotswdown ) ) THEN
          gotswdown = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'DSWF'
      END IF

      IF ( ( name == 'TERRAIN  ' ) .AND. ( .NOT. gotterrain ) ) THEN
          gotterrain = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'SHGT'
      END IF

      IF ( ( name == 'PBL HGT  ' ) .AND. ( .NOT. gotpblhgt ) ) THEN
          gotpblhgt = .TRUE.
          n0 = n0 + 1
          vchar0 (n0) = 'MXHT'
      END IF

      IF ( ( name == 'SIGMAH   ' ) .AND. ( .NOT. gotsig ) ) THEN
          gotsig = .TRUE.
          halfsig = data(:,1,1,1)

!     Standard sigma levels are converted to pseudo sigma levels where it is
!     assumed that ptop=0 to be consistent with old hysplit dispersion code. 
!     The updated code uses the actual pressure level data at each grid point.
          arlsig  = halfsig(:) + ptop * ( 1.0 - halfsig(:) ) / 1013.0

      END IF

      DEALLOCATE ( data )

    ELSE IF ( iflag == 2 ) THEN

      REWIND     (iutmm)
      DEALLOCATE ( bhi  )
      DEALLOCATE ( bhic )
      DEALLOCATE ( bhr  )
      DEALLOCATE ( bhrc )
      RETURN                

    ELSE

      PRINT*, '*** UNKNOWN MM5v3 DATA FLAG...IFLAG = ', iflag
      STOP

    ENDIF

  ENDDO v3data

!-------------------------------------------------------------------------------
! Error-handling section.
!-------------------------------------------------------------------------------

 8000 WRITE (6,9000) iutmm, TRIM(var), istat
      CALL abort
      STOP

 8100 WRITE (6,9100) iutmm, iflag, istat
      CALL abort
      STOP

 8200 WRITE (6,9200) iutmm, istat
      CALL abort
      STOP

 8300 WRITE (6,9300) iutmm, TRIM(var), istat
      CALL abort
      STOP

 9000 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** PROGRAM:  setvarb',                     &
              /, 1x, '***   ERROR READING FILE, UNIT = ', i3,     &
              /, 1x, '***   VARIABLE = ', a,                      &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9100 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** PROGRAM:  setvarb',                     &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   IFLAG = ', i3,                        &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9200 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** PROGRAM:  setvarb',                     &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   VARIABLE = SMALL HEADER',             &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

 9300 FORMAT (/, 1x, 70('*'),                                     &
              /, 1x, '*** PROGRAM:  setvarb',                     &
              /, 1x, '***   UNEXPECTED END-OF-FILE, UNIT = ', i3, &
              /, 1x, '***   VARIABLE = ', a,                      &
              /, 1x, '***   IOSTAT = ', i4,                       &
              /, 1x, 70('*'))

END SUBROUTINE setvarb
