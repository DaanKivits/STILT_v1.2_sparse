	SUBROUTINE CW2GLL (STCPRM, XLAT,XLONG, UE,VN, UG,VG)
!*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      DOUBLE PRECISION XPOLG,YPOLG,ALONG,SLONG,CLONG,ROT
      REAL STCPRM(9)
      ALONG = CSPANF( XLONG - STCPRM(2), -180., 180.)
      IF (XLAT.GT.89.) THEN
!*  NORTH POLAR METEOROLOGICAL ORIENTATION: "NORTH" ALONG PRIME MERIDIAN
        ROT = - STCPRM(1) * ALONG + XLONG - 180.
      ELSEIF (XLAT.LT.-89.) THEN
!*  SOUTH POLAR METEOROLOGICAL ORIENTATION: "NORTH" ALONG PRIME MERIDIAN
        ROT = - STCPRM(1) * ALONG - XLONG
      ELSE
        ROT = - STCPRM(1) * ALONG
      ENDIF
      SLONG = SIN( RADPDG * ROT )
      CLONG = COS( RADPDG * ROT )
      XPOLG = SLONG * STCPRM(5) + CLONG * STCPRM(6)
      YPOLG = CLONG * STCPRM(5) - SLONG * STCPRM(6)
!*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
      TEMP = YPOLG * UE + XPOLG * VN
      VG = YPOLG * VN - XPOLG * UE
      UG = TEMP
!* PERMITTING ROTATE OF WINDS IN PLACE
      RETURN
      END
