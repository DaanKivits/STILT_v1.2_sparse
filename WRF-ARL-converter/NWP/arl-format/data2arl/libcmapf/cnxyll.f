      SUBROUTINE CNXYLL (STCPRM, XI,ETA, XLAT,XLONG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
C  MAIN TRANSFORMATION ROUTINE FROM CANONICAL (EQUATOR-CENTERED,
C  RADIAN UNIT) COORDINATES
      include 'cmapf.fi'
      real stcprm(k_maparam)
      DOUBLE PRECISION GAMMA,TEMP,ARG1,ARG2,YMERC,ALONG,GXI,CGETA
      GAMMA = STCPRM(k_gama)
      CGETA = 1.D0 - GAMMA * ETA
      GXI = GAMMA * XI
C  CALCULATE EQUIVALENT MERCATOR COORDINATE
      ARG2 = ETA + (ETA * CGETA - GXI * XI)
      ARG1 = GAMMA * ARG2
      IF (ARG1 .GE. 1.0) THEN
C  DISTANCE TO NORTH (OR SOUTH) POLE IS ZERO (OR IMAGINARY ;) )
        XLAT = SIGN(90.,STCPRM(k_gama))
C        XLONG = STCPRM(k_reflon)
        XLONG = 90. + XLAT
C Change made 02/12/02 to acommodate WMO reporting conventions.  North
C pole is longitude 180., so "North" points to the Greenwich Meridian,
C South Pole is longitude 0. so "North" again points to the Greenwich
C Meridian.
        RETURN
      ENDIF
      IF (ABS(ARG1) .LT. .01) THEN
C  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
C  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
        TEMP = (ARG1 / (2. - ARG1) )**2
        YMERC = ARG2 / (2. - ARG1) * (1.    + TEMP *
     C                               (1./3. + TEMP *
     C                               (1./5. + TEMP *
     C                               (1./7. ))))
      ELSE
C CODE FOR MODERATE VALUES OF GAMMA
        YMERC = - LOG ( 1. - ARG1 ) /2. / GAMMA
      ENDIF
C  CONVERT YMERC TO LATITUDE
      TEMP = EXP( - ABS(YMERC) )
      XLAT = SIGN(ATAN2((1. - TEMP) * (1. + TEMP), 2. * TEMP), YMERC)
C  FIND LONGITUDES
      IF ( ABS(GXI) .LT. .01*CGETA ) THEN
C  CODE FOR GAMMA SMALL OR ZERO.  THIS AVOIDS ROUND-OFF ERROR OR DIVIDE-
C  BY ZERO IN THE CASE OF MERCATOR OR NEAR-MERCATOR PROJECTIONS.
        TEMP = ( GXI /CGETA )**2
        ALONG = XI / CGETA * (1.    - TEMP *
     C                         (1./3. - TEMP *
     C                         (1./5. - TEMP *
     C                         (1./7.   ))))
      ELSE
C CODE FOR MODERATE VALUES OF GAMMA
        ALONG = ATAN2( GXI, CGETA) / GAMMA
      ENDIF
      XLONG = SNGL(STCPRM(k_reflon) + DGPRAD * ALONG)
      XLAT = XLAT * DGPRAD
      RETURN
      END
