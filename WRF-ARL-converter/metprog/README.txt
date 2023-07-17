METEOROLOGY CONVERSION PROGRAMS
Last Updated: 29 January 2015


This directory contains the source code to compile the
five most popular programs to convert meteorological 
model output files into a format that can be read by
HYSPLIT. The three programs handle data from WRFARW
(requires NetCDF), GRIB2 formatted data (requires the
ECMWF grib_api library), and GRIB1 formatted data.
External librarys must be downloaded and installed prior
to compiling the programs. The library routines included
here are required to create the HYSPLIT formatted
meteorological input data files.  Each directory contains
a Makefile and further instructions. Conversion programs
for MM5 and RAMS are provided for consistency, however
these programs are no longer used nor have they been
tested with the current HYSPLIT library.

_________________________________________________________
NETCDF Conversion - arw2arl

	Converts WRF-ARW NetCDF output files.
	Requires the installation of the NetCDF library.


_________________________________________________________
GRIB1 Conversion - grib2arl
 
	Converts GRIB version 1 data files.
	Input data must be on a latitude-longitude grid.
	Required subroutines are in the library directory.
	
_________________________________________________________
GRIB2 Conversion - api2arl

	Converts GRIB version 2 data files. 
	Input data must be on a latitude-longitude grid.
	Requres installation of the ECMWF grib_api libraries. 
	http://www.ecmwf.int/products/data/software/download/grib_api.html

	Versions 

	v1 : original version processes pressure level data
             from the NCEP NOMADS server

	v2 : modified to handle the hybrid sigma coordinate
             on a Lambert Conformal projection and where the
             vertical index increases with height

	v3 : revised to handle 0.5 degree data archived by
             from NOAA (pressure) and ECMWF (hybrid) and
             where the vertical index decreases with height

	v4 : revised to handle GSD HRRR data on pressure
             level surfaces; uncomment sections for precipitation
             accumulations at intervals different than one hour
