
      SUBROUTINE CCRVLL (STRCMP, XLAT,XLONG, GX,GY)
!*  WRITTEN ON 9/20/94 BY Dr. Albion Taylor  NOAA / OAR / ARL
      PARAMETER (REARTH=6371.2)
      PARAMETER (PI=3.14159265358979,RADPDG=PI/180,DGPRAD=180/PI)
      DOUBLE PRECISION XPOLG,YPOLG,TEMP,ALONG,SLONG,CLONG,CTEMP
      REAL STRCMP(9)
      ALONG = CSPANF( XLONG - STRCMP(2), -180., 180.)
      SLONG = SIN( RADPDG * STRCMP(1) * ALONG)
      CLONG = COS( RADPDG * STRCMP(1) * ALONG)
      XPOLG = - SLONG * STRCMP(5) + CLONG * STRCMP(6)
      YPOLG = CLONG * STRCMP(5) + SLONG * STRCMP(6)
        TEMP = SIN(RADPDG * XLAT)
        CTEMP = COS(RADPDG * XLAT)
        CURV = (STRCMP(1) - TEMP) / CTEMP / REARTH
        GX = CURV * XPOLG
        GY = CURV * YPOLG
      RETURN
      END
