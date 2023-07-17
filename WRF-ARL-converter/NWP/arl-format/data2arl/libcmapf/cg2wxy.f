      SUBROUTINE CG2WXY (STCPRM, X,Y, UG,VG, UE,VN)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,XI0,ETA0
      XI0 = ( X - STCPRM(k_x0) ) * STCPRM(k_gdszeq) / REARTH
      ETA0 = ( Y - STCPRM(k_y0) ) * STCPRM(k_gdszeq) /REARTH
      XI = XI0 * STCPRM(k_crot) - ETA0 * STCPRM(k_srot)
      ETA = ETA0 * STCPRM(k_crot) + XI0 * STCPRM(k_srot)
      RADIAL = 2. * ETA - STCPRM(k_gama) * (XI*XI + ETA*ETA)
      IF (RADIAL.GT.STCPRM(k_n1dgr)) THEN
C*  CASE NORTH OF 89 DEGREES.  METEOROLOGICAL WIND DIRECTION DEFINITION
C*      CHANGES.
         CALL CNXYLL(STCPRM, XI,ETA, XLAT,XLONG)
         CALL CG2WLL(STCPRM, XLAT,XLONG, UG,VG, UE,VN)
      ELSE IF (RADIAL.LT.STCPRM(k_s1dgr)) THEN
C*  CASE SOUTH OF -89 DEGREES.  METEOROLOGICAL WIND DIRECTION DEFINITION
C*      CHANGES.
         CALL CNXYLL(STCPRM, XI,ETA, XLAT,XLONG)
         CALL CG2WLL(STCPRM, XLAT,XLONG, UG,VG, UE,VN)
      ELSE
C* NORMAL CASE.  METEOROLOGICAL DIRECTION OF WIND RELATED TO TRUE NORTH.
         XPOLG = STCPRM(k_srot) - STCPRM(k_gama) * XI0
         YPOLG = STCPRM(k_crot) - STCPRM(k_gama) * ETA0
         TEMP = SQRT ( XPOLG ** 2 + YPOLG ** 2 )
         XPOLG = XPOLG / TEMP
         YPOLG = YPOLG / TEMP
C*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
         TEMP = YPOLG * UG - XPOLG * VG
         VN = YPOLG * VG + XPOLG * UG
         UE = TEMP
C* PERMITTING ROTATE OF WINDS IN PLACE
      END IF
      RETURN
      END
