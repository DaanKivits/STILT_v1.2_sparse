      REAL FUNCTION CGSZXY (STCPRM, X,Y)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      real stcprm(k_maparam)
      PARAMETER (ALMST1=.99999)
      DOUBLE PRECISION YMERC,EFACT
      DOUBLE PRECISION XI0,ETA0,XI,ETA
      XI0 = ( X - STCPRM(k_x0) ) * STCPRM(k_gdszeq) / REARTH
      ETA0 = ( Y - STCPRM(k_y0) ) * STCPRM(k_gdszeq) /REARTH
      XI = XI0 * STCPRM(k_crot) - ETA0 * STCPRM(k_srot)
      ETA = ETA0 * STCPRM(k_crot) + XI0 * STCPRM(k_srot)
      RADIAL = 2. * ETA - STCPRM(k_gama) * (XI*XI + ETA*ETA)
      EFACT = STCPRM(k_gama) * RADIAL
      IF (EFACT .GT. ALMST1) THEN
         IF (STCPRM(k_gama).GT.ALMST1) THEN
	         CGSZXY = 2. * STCPRM(k_gdszeq)
         ELSE
            CGSZXY = 0.
         ENDIF
         RETURN
      ENDIF
      IF (ABS(EFACT) .LT. 1.E-2) THEN
         TEMP = (EFACT / (2. - EFACT) )**2
         YMERC = RADIAL / (2. - EFACT) * (1.    + TEMP *
     C				          (1./3. + TEMP *
     C				          (1./5. + TEMP *
     C				          (1./7. ))))
      ELSE
         YMERC = - LOG( 1. - EFACT ) /2. /STCPRM(k_gama)
      ENDIF
      IF (YMERC .GT. 6.) THEN
         IF (STCPRM(k_gama) .GT. ALMST1) THEN
            CGSZXY = 2. * STCPRM(k_gdszeq)
         ELSE
            CGSZXY = 0.
         ENDIF
      ELSE IF (YMERC .LT. -6.) THEN
         IF (STCPRM(k_gama) .LT. -ALMST1) THEN
            CGSZXY = 2. * STCPRM(k_gdszeq)
         ELSE
            CGSZXY = 0.
         ENDIF
      ELSE
         EFACT = EXP(YMERC)
         CGSZXY = 2. * STCPRM(k_gdszeq) * EXP (STCPRM(k_gama) * YMERC)
     C				 / (EFACT + 1./EFACT)
      ENDIF
      RETURN
      END
