			NCEP/NCAR Reanalysis Data Archive 
				Revised: 04 Apr 2001

________________________________________________________________________________

Overview -

The NCEP/NCAR Reanalysis Project is a joint project between the National 
Centers for Environmental Prediction (NCEP, formerly "NMC") and the National 
Center for Atmospheric Research (NCAR). The goal of this joint effort is to 
produce new atmospheric analyses using historical data (1948 onwards) and as 
well to produce analyses of the current atmospheric state (Climate Data 
Assimilation System, CDAS).  Until recently, the meteorological community has 
had to use analyses that supported the real-time weather forecasting. These 
analyses are very inhomogeneous in time as there have been big improvements in 
the data assimilation systems. The quality and utility of the re-analyses 
should be superior to NCEP's original analyses because:

* a state-of-the-art data assimilation is used 
* more observations are used 
* quality control has been improved 
* the model/data assimilation procedure remains unchanged during the project 
* many more fields are being saved 
* global (some older analyses were hemispheric) 
* better vertical resolution (stratosphere) 

More information about the reanalysis project and data are available from
several sources:

http://wesley.wwb.noaa.gov/reanalysis.html
http://www.cdc.noaa.gov/cdc/data.nmc.reanalysis.html

________________________________________________________________________________

Availability -

A subset of this data is available from ARL in a format suitable for transport 
and dispersion calculations using HYSPLIT through READY by selecting 
"Reanalysis" in the meteorological data set selection pull-down menu.  Data 
files are identified by the following syntax:

R{S|P}{YEAR}{MONTH}.{gbl|usa|tbd}

Where R indicates "Reanalysis", S or P indicates that the data are on Sigma or 
Pressure surfaces, YEAR is a four digit year, and MONTH is a two digit month. 
The file suffix identifies the projection as either the 2.5 degree global 
latitude-longitude projection (gbl), or a regional conformal map projection 
(such as over the "usa").  Other regional projections are "to be determined 
(tbd)" later.  The projection details are encoded in the file's index record 
and are processed by HYSPLIT during trajectory or dispersion computations.

The sigma level data were obtained from NCEP's internal spectral coefficient
archive. The pressure level data were obtained from the NOAA-CIRES Climate 
Diagnostics Center, Boulder, Colorado, USA, through their Web site at 
http://www.cdc.noaa.gov/. 

Due to storage limitations only a fraction of the 40 year reanalyis period is 
available for on-line computations. Other product distribution mechanisms are
under investigation.  Until more automated procedures are in-place do not ask 
for additional data to be placed on-line. All years of Reanalysis data (netCDF
format) may be obtained on-line through the CDC web site (see CDC link above)
as yearly files per variable.  Eleven files are required to create the HYSPLIT
compatible input data set. The required file listing for 1987 is shown below:

	512 blks  Data set name 

	521830088 air.1987.nc
	 30709096 air.sig995.1987.nc
	521830096 hgt.1987.nc
	368354864 omega.1987.nc
	 52714616 prate.sfc.gauss.1987.nc
	 30709140 pres.sfc.1987.nc
	245574696 rhum.1987.nc
	521829984 uwnd.1987.nc
	 30709080 uwnd.sig995.1987.nc
	521829984 vwnd.1987.nc
	 30709080 vwnd.sig995.1987.nc

A conversion program from CDC netCDF format to ARL HYSPLIT compatible format is
available upon requrest.  The compilation procedures are rather complex and 
the processing requires substantial storage space.  No support can be provided.

________________________________________________________________________________

Additional Data Set Details

  Pressure Level Data

	2.5 degree latitude-longitude global grid
	144x73 points from 90N-90S, 0E-357.5E 
	1/1/1948 - present with output every 6 hours 
	Levels (hPa): 1000,925,850,700,600,500,400,300,250,200,
                      150,100,70,50,30,20,10 
	Surface or near the surface (.995 sigma level) winds and temperature 
	Precipitation 

	Model Type:            LAT-LON
	Vert Coord:            2
	Numb X pt:           144
	Numb Y pt:            73
	Numb Levels:          18
	Sfc Variables:         5 PRSS T02M U10M V10M TPP6
	Upper Levels:          6 HGTS TEMP UWND VWND WWND RELH


  Sigma Level Data

	The spectral coefficients on 28 model sigma surfaces were processed
	to obtain required fields 4 per day on a global Gaussian grid of 
	1.875 degree resolution.  A regional sub-grid covering the continental 
	US and Canada was extracted.

	Current USA spatial domain: 21.9N 127.5W to 60.0N 52.5W
	Output every 6 hours
	Levels: .995,.982,.964,.943,.916,.884,.846,.801,.751,.694,.633,
                .568,.502,.436,.372,.312,.258,.210,.168,.133,.103,.078,
                .058,.042,.029,.018,.010,.003 

	MODEL TYPE:         MERCATOR
	VERT COORD:         1
	POLE LAT:           21.904
	POLE LON:           -127.5
	REF LAT:            21.904
	REF LON:            -127.5
	REF GRID:           136.5
	ORIENTATION:        0.
	CONE ANGLE:         0.
	SYNC X:             1.
	SYNC Y:             1.
	SYNC LAT:           21.904
	SYNC LON:           -127.5
	NUMB X:             57
	NUMB Y:             41
	NUMB LEVELS:        29
	SFC VARIABLES:      01 PRSS
	UPPER LEVELS:       05 TEMP SPHU UWND VWND WWND

