MODULE consts

  REAL,    PARAMETER :: cp     =  1004.0    ! J/kg/K; specific heat
  REAL,    PARAMETER :: dummy  =   100.0    ! dummy value
  REAL,    PARAMETER :: g      =     9.81   ! m/s**2; gravity
  REAL,    PARAMETER :: pi     =     3.1415926535
  REAL,    PARAMETER :: r      =   287.04   ! J/kg/K; gas constant
  REAL,    PARAMETER :: t0     =   273.15   ! 0 Celsius in Kelvin
  REAL,    PARAMETER :: vkar   =     0.4    ! Von Karman constant
  REAL,    PARAMETER :: xmiss  =  -999.99   ! missing data
  REAL,    PARAMETER :: zsfct  =    10.0    ! m; shelter height for T
  REAL,    PARAMETER :: zsfcv  =    10.0    ! m; tower height for wind

  REAL,    PARAMETER :: degran = 180.0 / pi ! degrees to radians
  REAL,    PARAMETER :: rovcp  = r / cp

END MODULE consts
