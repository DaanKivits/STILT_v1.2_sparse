      SUBROUTINE CC2GLL (STCPRM, XLAT,XLONG, UE,VN, UG,VG)
C*  WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      DOUBLE PRECISION XPOLG,YPOLG,ALONG,SLONG,CLONG,ROT
      ALONG = CSPANF( XLONG - STCPRM(k_reflon), -180., 180.)
C* Revised 2/12/02 to allow cartographic wind vector transformations everywhere
C* with rotation to nominal longitudes at the poles, to match U,V values on a
C* Lat-Lon grid.
      ROT = - STCPRM(k_gama) * ALONG
      SLONG = SIN( RADPDG * ROT )
      CLONG = COS( RADPDG * ROT )
      XPOLG = SLONG * STCPRM(k_crot) + CLONG * STCPRM(k_srot)
      YPOLG = CLONG * STCPRM(k_crot) - SLONG * STCPRM(k_srot)
C*  CHANGE MADE 3/9/99 TO ALLOW UG,VG TO HAVE SAME STORAGE AS UE,VN
      TEMP = YPOLG * UE + XPOLG * VN
      VG = YPOLG * VN - XPOLG * UE
      UG = TEMP
C* PERMITTING ROTATE OF WINDS IN PLACE
      RETURN
      END
