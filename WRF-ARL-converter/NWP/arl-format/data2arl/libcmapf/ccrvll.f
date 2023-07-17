      SUBROUTINE CCRVLL (STCPRM, XLAT,XLONG, GX,GY)
C*  WRITTEN ON 9/20/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,ALONG,SLONG,CLONG,CTEMP
      ALONG = CSPANF( XLONG - STCPRM(k_reflon), -180., 180.)
      SLONG = SIN( RADPDG * STCPRM(k_gama) * ALONG)
      CLONG = COS( RADPDG * STCPRM(k_gama) * ALONG)
      XPOLG = - SLONG * STCPRM(k_crot) + CLONG * STCPRM(k_srot)
      YPOLG = CLONG * STCPRM(k_crot) + SLONG * STCPRM(k_srot)
      TEMP = SIN(RADPDG * XLAT)
      CTEMP = COS(RADPDG * XLAT)
      CURV = (STCPRM(k_gama) - TEMP) / CTEMP / REARTH
      GX = CURV * XPOLG
      GY = CURV * YPOLG
      RETURN
      END
