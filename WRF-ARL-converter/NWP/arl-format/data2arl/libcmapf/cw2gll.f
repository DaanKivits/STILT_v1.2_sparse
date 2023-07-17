      SUBROUTINE CW2GLL (STCPRM, XLAT,XLONG, UE,VN, UG,VG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      ALONG = CSPANF( XLONG - STCPRM(k_reflon), -180., 180.)
      IF (XLAT.GT.89.) THEN
C*  NORTH POLAR METEOROLOGICAL ORIENTATION: "NORTH" ALONG PRIME MERIDIAN
         ROT = - STCPRM(k_gama) * ALONG + XLONG - 180.
      ELSEIF (XLAT.LT.-89.) THEN
C*  SOUTH POLAR METEOROLOGICAL ORIENTATION: "NORTH" ALONG PRIME MERIDIAN
         ROT = - STCPRM(k_gama) * ALONG - XLONG
      ELSE
         ROT = - STCPRM(k_gama) * ALONG
      ENDIF
      SLONG = SIN( RADPDG * ROT )
      CLONG = COS( RADPDG * ROT )
      XPOLG = SLONG * STCPRM(k_crot) + CLONG * STCPRM(k_srot)
      YPOLG = CLONG * STCPRM(k_crot) - SLONG * STCPRM(k_srot)
C*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
      TEMP = YPOLG * UE + XPOLG * VN
      VG = YPOLG * VN - XPOLG * UE
      UG = TEMP
C* PERMITTING ROTATE OF WINDS IN PLACE
      RETURN
      END
