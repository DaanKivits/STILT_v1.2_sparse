      SUBROUTINE CC2GXY (STCPRM, X,Y, UE,VN, UG,VG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,XI0,ETA0
      XI0 = ( X - STCPRM(k_x0) ) * STCPRM(k_gdszeq) / REARTH
      ETA0 = ( Y - STCPRM(k_y0) ) * STCPRM(k_gdszeq) /REARTH
      XPOLG = STCPRM(k_srot) - STCPRM(k_gama) * XI0
      YPOLG = STCPRM(k_crot) - STCPRM(k_gama) * ETA0
      TEMP = SQRT ( XPOLG ** 2 + YPOLG ** 2 )
C* Revised 2/12/02 to allow cartographic wind vector transformations everywhere
C* except at the poles, with WMO conventions only at the poles.
      IF (TEMP .LE. 0.1e-3) THEN
         ANG = RADPDG * STCPRM(k_gama)* STCPRM(k_reflon)
         IF (STCPRM(k_gama) .GT. 0.) THEN
C* North Pole case; vector directed along the Greenwich meridian
	         XPOLG = - SIN(ANG)
            YPOLG = - COS(ANG)
         ELSE
C* South Pole case; vector directed along the Greenwich meridian
            XPOLG = SIN(ANG)
            YPOLG = COS(ANG)
         ENDIF
      ELSE
         XPOLG = XPOLG / TEMP
         YPOLG = YPOLG / TEMP
      ENDIF
C*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
      TEMP = YPOLG * UE + XPOLG * VN
      VG = YPOLG * VN - XPOLG * UE
      UG = TEMP
C* PERMITTING ROTATE OF WINDS IN PLACE
      RETURN
      END
