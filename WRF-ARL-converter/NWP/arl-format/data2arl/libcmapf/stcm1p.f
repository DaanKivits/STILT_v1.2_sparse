      SUBROUTINE STCM1P(STCPRM, X1,Y1, XLAT1,XLONG1,
     C XLATG,XLONGG, GRIDSZ, ORIENT)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      stcprm(k_x0) = 0.
      stcprm(k_y0) = 0.
      TURN = RADPDG * (ORIENT - STCPRM(k_gama) *
     C            CSPANF(XLONGG - STCPRM(k_reflon), -180., 180.) )
      STCPRM (k_crot) = COS (TURN)
      STCPRM (k_srot) = - SIN (TURN)
      STCPRM (k_gdszeq) = 1.
      STCPRM (k_gdszeq) = GRIDSZ * STCPRM(k_gdszeq)
     C             / CGSZLL(STCPRM, XLATG, STCPRM(k_reflon))
	CALL CLL2XY (STCPRM, XLAT1,XLONG1, X1A,Y1A)
      STCPRM(k_x0) = STCPRM(k_x0) + X1 - X1A
      STCPRM(k_y0) = STCPRM(k_y0) + Y1 - Y1A
      RETURN
      END
