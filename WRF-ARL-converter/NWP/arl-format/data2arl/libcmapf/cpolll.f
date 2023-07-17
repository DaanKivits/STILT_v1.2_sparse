      SUBROUTINE CPOLLL (STCPRM, XLAT,XLONG, ENX,ENY,ENZ)
C*  WRITTEN ON 11/23/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,ALONG,SLONG,CLONG,ROT
      ALONG = CSPANF( XLONG - STCPRM(k_reflon), -180., 180.)
      ROT = - STCPRM(k_gama) * ALONG
      SLONG = SIN( RADPDG * ROT )
      CLONG = COS( RADPDG * ROT )
      XPOLG = SLONG * STCPRM(k_crot) + CLONG * STCPRM(k_srot)
      YPOLG = CLONG * STCPRM(k_crot) - SLONG * STCPRM(k_srot)
      CLAT = COS(RADPDG * XLAT)
      ENX = CLAT * XPOLG
      ENY = CLAT * YPOLG
      ENZ = SIN(RADPDG * XLAT)
      RETURN
      END
