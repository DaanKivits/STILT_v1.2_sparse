      SUBROUTINE CPOLXY (STCPRM, X,Y, ENX,ENY,ENZ)
C*  WRITTEN ON 11/26/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOL,YPOL,TEMP,XI0,ETA0,XI,ETA,RADIAL
      DOUBLE PRECISION TEMP2,YMERC,ARG,OARG,CLAT
      XI0 = ( X - STCPRM(k_x0) ) * STCPRM(k_gdszeq) / REARTH
      ETA0 = ( Y - STCPRM(k_y0) ) * STCPRM(k_gdszeq) /REARTH
      XI = XI0 * STCPRM(k_crot) - ETA0 * STCPRM(k_srot)
      ETA = ETA0 * STCPRM(k_crot) + XI0 * STCPRM(k_srot)
      RADIAL = 2. * ETA -  STCPRM(k_gama) * (XI*XI + ETA*ETA)
      TEMP = STCPRM(k_gama) * RADIAL
      IF (TEMP .GE. 1.) THEN
         ENX = 0.
         ENY = 0.
         ENZ = SIGN(1.,STCPRM(k_gama))
         RETURN
      ENDIF
      IF (ABS(TEMP).LT.1.E-2) THEN
         TEMP2 = (TEMP / (2. - TEMP))**2
         YMERC = RADIAL / (2. - TEMP) * (1. + TEMP2 *
     C					 (1./3. + TEMP2 *
     C					 (1./5. + TEMP2 *
     C					 (1./7.))))
      ELSE
         YMERC = -.5 * LOG(1. - TEMP) / STCPRM(k_gama)
      ENDIF
      ARG = EXP( YMERC )
      OARG = 1./ARG
      CLAT = 2./(ARG + OARG)
      ENZ = (ARG - OARG) * CLAT /2.
      TEMP = CLAT / SQRT(1. - TEMP)
      XPOL = - XI * STCPRM(k_gama) * TEMP
      YPOL = (1. - ETA * STCPRM(k_gama) ) * TEMP
      ENX = XPOL * STCPRM(k_crot) + YPOL * STCPRM(k_srot)
      ENY = YPOL * STCPRM(k_crot) - XPOL * STCPRM(k_srot)
      RETURN
      END
