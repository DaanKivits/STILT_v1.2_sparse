      SUBROUTINE CXY2LL (STCPRM, X,Y, XLAT,XLONG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      real stcprm(k_maparam)
      XI0 = ( X - STCPRM(k_x0) ) * STCPRM(k_gdszeq) / REARTH
      ETA0 = ( Y - STCPRM(k_y0) ) * STCPRM(k_gdszeq) /REARTH
      XI = XI0 * STCPRM(k_crot) - ETA0 * STCPRM(k_srot)
      ETA = ETA0 * STCPRM(k_crot) + XI0 * STCPRM(k_srot)
      CALL CNXYLL(STCPRM, XI,ETA, XLAT,XLONG)
      XLONG = CSPANF(XLONG, -180., 180.)
      RETURN
      END
