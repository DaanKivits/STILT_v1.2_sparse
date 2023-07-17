      SUBROUTINE STCM2P(STCPRM, X1,Y1, XLAT1,XLONG1,
     C X2,Y2, XLAT2,XLONG2)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      real stcprm(k_maparam)
      stcprm(k_x0) = 0.
      stcprm(k_y0) = 0.
      stcprm(k_crot) = 1.
      stcprm(k_srot) = 0.
      STCPRM (k_gdszeq) = 1.
      CALL CLL2XY (STCPRM, XLAT1,XLONG1, X1A,Y1A)
      CALL CLL2XY (STCPRM, XLAT2,XLONG2, X2A,Y2A)
      DEN = SQRT( (X1 - X2)**2 + (Y1 - Y2)**2 )
      DENA = SQRT( (X1A - X2A)**2 + (Y1A - Y2A)**2 )
      STCPRM(k_crot) = ((X1A - X2A)*(X1 - X2) + (Y1A - Y2A)*(Y1 - Y2))
     C  /DEN /DENA
      STCPRM(k_srot) = ((Y1A - Y2A)*(X1 - X2) - (X1A - X2A) * (Y1 - Y2))
     C  /DEN /DENA
      STCPRM (k_gdszeq) = STCPRM(k_gdszeq) * DENA / DEN
      CALL CLL2XY (STCPRM, XLAT1,XLONG1, X1A,Y1A)
      STCPRM(k_x0) = STCPRM(k_x0) + X1 - X1A
      STCPRM(k_y0) = STCPRM(k_y0) + Y1 - Y1A
      RETURN
      END
