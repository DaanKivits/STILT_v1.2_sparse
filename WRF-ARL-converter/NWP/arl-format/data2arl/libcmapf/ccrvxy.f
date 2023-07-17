      SUBROUTINE CCRVXY (STCPRM, X,Y, GX,GY)
C*  WRITTEN ON 9/20/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,YMERC,EFACT,CURV
      TEMP = STCPRM(k_gama) * STCPRM(k_gdszeq) /REARTH
      XPOLG = STCPRM(k_srot) + TEMP * (STCPRM(k_x0) - X)
      YPOLG = STCPRM(k_crot) + TEMP * (STCPRM(k_y0) - Y)
      TEMP = SQRT ( XPOLG ** 2 + YPOLG ** 2 )
      IF (TEMP.GT.0.) THEN
         YMERC = - LOG( TEMP) /STCPRM(k_gama)
         EFACT = EXP(YMERC)
         CURV = ( (STCPRM(k_gama) - 1.D0) * EFACT +
     A            (STCPRM(k_gama) + 1.D0) / EFACT )
     B           * .5D0 / REARTH
         GX = XPOLG * CURV / TEMP
         GY = YPOLG * CURV / TEMP
      ELSE
         IF (ABS(STCPRM(k_gama)) .EQ. 1.) THEN
            GX = 0.
            GY = 0.
         ELSE
            GX = 1./REARTH
            GY = 1./REARTH
         ENDIF
      ENDIF
      RETURN
      END
