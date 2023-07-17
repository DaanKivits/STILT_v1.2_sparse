SUBROUTINE setupv3 (fname, ix, jx, kx, ptop, p0s, ts0, tlp, tiso,    &
                    nicg, njcg, proj, gnum, dist, clat, clon,        &
                    tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2) 

!-------------------------------------------------------------------------------
! Name:     Set Up the Domain Attributes from MM5 Version 3 Output
! Purpose:  Establishes bounds for post-processing.
! Method:   Extract IX, JX, and KX and rewind unit.
! Revised:  08 Mar 2000  Original version for MM5v3.  (TLO)
!           09 Aug 2000  mm5toarl application.  (RRD) 
!           09 Mar 2004  pass through domain number (RRD)
!-------------------------------------------------------------------------------

  USE file

  IMPLICIT NONE

  INTEGER,       ALLOCATABLE   :: bhi       ( : , : )
  REAL,          ALLOCATABLE   :: bhr       ( : , : )
  INTEGER                      :: iflag
  INTEGER                      :: index
  INTEGER                      :: istat
  CHARACTER*80,  INTENT(IN)    :: fname
  INTEGER,       INTENT(OUT)   :: ix
  INTEGER,       INTENT(OUT)   :: jx
  INTEGER,       INTENT(OUT)   :: kx
  INTEGER,       PARAMETER     :: numprogs  = 20
  INTEGER,       PARAMETER     :: numvalsi  = 50
  INTEGER,       PARAMETER     :: numvalsr  = 20
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

!-------------------------------------------------------------------------------
! Allocate necessary arrays.
!-------------------------------------------------------------------------------

  ALLOCATE ( bhi (numvalsi, numprogs) )
  ALLOCATE ( bhr (numvalsr, numprogs) )

  WRITE (6,6000)

!-------------------------------------------------------------------------------
! Read MM5 header for this domain.  Note that the real header information (BHR)
! and character header information (BHIC and BHRC) are not read or used here.
!-------------------------------------------------------------------------------

! Required when transferring data between some LINUX and UNIX platforms
! OPEN (iutmm, FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
  OPEN (iutmm, FILE=FNAME,FORM='UNFORMATTED')

  READ (iutmm, IOSTAT=istat, ERR=8000, END=8100) iflag
  
  IF ( iflag == 0 ) THEN   
    READ (iutmm, IOSTAT=istat, ERR=8000, END=8100) bhi, bhr
  ELSE
    GO TO 8200
  ENDIF

!-------------------------------------------------------------------------------
! Extract IX, JX, and KX
!-------------------------------------------------------------------------------

  index = bhi(1,1)
  IF ( index /= 11 ) GOTO 8300  ! output must be from mm5 v3

  ix = bhi(16,1)            ! domain grid dimension
  jx = bhi(17,1)            ! for data grid
  kx = bhi(12,index)        ! number of layers in output

  WRITE (6,6100) ix, jx, kx

!-------------------------------------------------------------------------------
! Extract some additional pressure/height scaling variables from header
!-------------------------------------------------------------------------------

  ptop = bhr(2,2) / 100.0   ! convert from Pa to hPa (mb)
  p0s  = bhr(2,5) / 1000.0  ! convert from Pa to kPa
  ts0  = bhr(3,5)           ! base state sea level temperature
  tlp  = bhr(4,5)           ! base state lapse rate
  tiso = bhr(5,5)           ! base state isothermal stratospheric temperature

!-------------------------------------------------------------------------------
! Extract geographic grid reference information for coarse domain
!-------------------------------------------------------------------------------

  nicg = bhi(5,1)           ! number of i points on coarse grid (S->N) 
  njcg = bhi(6,1)           ! number of j points on coarse grid (W->E)
  proj = bhi(7,1)           ! projection type 1=lambert 2=polar 3=mercator
  gnum = bhi(13,1)*10 &
        +bhi(14,1)          ! domain and mother grid ids
  dist = bhr(1,1)/1000.0    ! grid distance in kilometers  
  clat = bhr(2,1)           ! grid center latitude
  clon = bhr(3,1)           ! grid center longitude
  tru1 = bhr(5,1)           ! true latitude 1
  tru2 = bhr(6,1)           ! true latitude 2
  pole = bhr(7,1)           ! pole position latitude

  cgi1 = bhr(10,1)          ! i location in coarse domain of 1,1
  cgj1 = bhr(11,1)          ! j location in coarse domain of 1,1
  cgi2 = bhr(12,1)          ! i location in coarse domain of ix,jx
  cgj2 = bhr(13,1)          ! j location in coarse domain of ix,jx

!-------------------------------------------------------------------------------
! Rewind model output file
!-------------------------------------------------------------------------------

  REWIND (iutmm)

!-------------------------------------------------------------------------------
! Deallocate arrays.
!-------------------------------------------------------------------------------

  DEALLOCATE ( bhi )
  DEALLOCATE ( bhr )

  RETURN

!-------------------------------------------------------------------------------
! Format statements.
!-------------------------------------------------------------------------------

 6000 FORMAT (/, 1x, '- SUBROUTINE SETUPV3 - READING HEADER')
 6100 FORMAT (3x, 'GRID DIMENSIONS (I,J,K) ', i4, 1x, i4, 1x, i3)

!-------------------------------------------------------------------------------
! Error-handling section.
!-------------------------------------------------------------------------------

 8000 WRITE (6,9000) iutmm, istat
      CALL abort
      RETURN

 8100 WRITE (6,9100) iutmm, istat
      CALL abort
      RETURN

 8200 WRITE (6,9200) iflag

 8300 WRITE (6,9300) index
      CALL abort
      RETURN

 9000 FORMAT (/, 1x, 70('*'),                             &
              /, 1x, '*** SUBROUTINE: SETUPV3',           &
              /, 1x, '***   ERROR READING HEADER RECORD', &
              /, 1x, '***   UNIT = ', i3,                 &
              /, 1x, '***   IOSTAT = ', i4,               &
              /, 1x, 70('*'))

 9100 FORMAT (/, 1x, 70('*'),                                              &
              /, 1x, '*** SUBROUTINE: SETUPV3',                            &
              /, 1x, '***   UNEXPECTED END-OF-FILE REACHED ON UNIT ', i3,  &
              /, 1x, '***   IOSTAT = ', i5,                                &
              /, 1x, '***   VERIFY THAT THE FILE EXISTS!!!'                &
              /, 1x, 70('*'))

 9200 FORMAT (/, 1x, 70('*'),                                              &
              /, 1x, '*** SUBROUTINE: SETUPV3',                            &
              /, 1x, '***   UNEXPECTED FLAG FOUND IN VERSION 3 HEADER',    &
              /, 1X, '***   IFLAG = ', i3,                                 &
              /, 1x, 70('*'))

 9300 FORMAT (/, 1x, 70('*'),                                    &
              /, 1x, '*** SUBROUTINE: SETUPV3'                   &
              /, 1x, '***   INAPPROPRIATE INPUT FILE'            &
              /, 1x, '***   MUST BE INTERPF OUTPUT (MMINPUT)'    &
              /, 1x, '***   POSSIBLE UNIX/PC BYTE-SWAP'          &
              /, 1x, '***   INDEX IS ', i2,                      &
              /, 1x, 70('*'))

END SUBROUTINE setupv3
