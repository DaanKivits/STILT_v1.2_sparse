      REAL FUNCTION CGSZLL (STCPRM, XLAT,XLONG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      real stcprm(k_maparam)
      DOUBLE PRECISION SLAT,YMERC,EFACT
      IF (XLAT .GT. 89.995) THEN
C* CLOSE TO NORTH POLE
         IF (STCPRM(k_gama) .GT. 0.9999) THEN
C* AND TO GAMMA == 1.
            CGSZLL = 2. * STCPRM(k_gdszeq)
            RETURN
         ENDIF
         EFACT = COS(RADPDG * XLAT)
         IF (EFACT .LE. 0.) THEN
            CGSZLL = 0.
            RETURN
         ELSE
            YMERC = - LOG( EFACT /(1. + SIN(RADPDG * XLAT)))
         ENDIF
      ELSE IF (XLAT .LT. -89.995) THEN
C* CLOSE TO SOUTH POLE
         IF (STCPRM(k_gama) .LT. -0.9999) THEN
C* AND TO GAMMA == -1.0
            CGSZLL = 2. * STCPRM(k_gdszeq)
            RETURN
         ENDIF
         EFACT = COS(RADPDG * XLAT)
         IF (EFACT .LE. 0.) THEN
            CGSZLL = 0.
            RETURN
         ELSE
            YMERC = LOG( EFACT /(1. - SIN(RADPDG * XLAT)))
         ENDIF
      ELSE
         SLAT = SIN(RADPDG * XLAT)
         YMERC = LOG((1. + SLAT) / (1. - SLAT))/2.
C	EFACT = EXP(YMERC)
C	CGSZLL = 2. * STCPRM(7) * EXP (STCPRM(1) * YMERC)
C     C			 / (EFACT + 1./EFACT)
      ENDIF
      CGSZLL = STCPRM(k_gdszeq) * COS(RADPDG * XLAT) *
     C          EXP(STCPRM(k_gama) *YMERC)
      RETURN
      END
