      SUBROUTINE CLL2XY (STCPRM, XLAT,XLONG, X,Y)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      real stcprm(k_maparam)
      CALL CNLLXY(STCPRM, XLAT,XLONG, XI,ETA)
      X = STCPRM(k_x0) + REARTH/STCPRM(k_gdszeq) *
     C		 (XI * STCPRM(k_crot) + ETA * STCPRM(k_srot) )
      Y = STCPRM(k_y0) + REARTH/STCPRM(k_gdszeq) *
     C		 (ETA * STCPRM(k_crot) - XI * STCPRM(k_srot) )
      RETURN
      END
