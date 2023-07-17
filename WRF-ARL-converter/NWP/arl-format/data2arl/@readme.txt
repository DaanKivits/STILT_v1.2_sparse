DIRECTORY: /data2arl                      Last Revised: 17 Nov 2003
___________________________________________________________________

This directory contains various directories with example programs to 
convert meteorological data in various formats to the format (ARL 
packed) that Hysplit can read. The use of the packing routines is 
described in the User's Guide.  Routines may require some customiztion
depending upon the structure of the input data files.  

Briefly a ARL packed data file consists of a series of fixed length 
records, one for each meteorological variable. The records are 
arranged in time series, with surface fields followed by upper air
fields.  Each record contains a 50 byte ascii header that provides
information about the date, time, variable, and packing constants.
This is followed by the packed binary data, one byte per grid point.
Each byte represents the difference in value between that point and
the previous grid point.  A group of records at the same time is 
preceeded by an index record that describes all the variables, levels,
and grid information for that time period.


The following programs are available: 

	afwa -		Air Force Weather Agency's MM5 GRIB output
			afwa2arl: AWIPS-212 grid 
			pNA45: 45 km Lambert Conformal projection
			pNA15: 15 km resolution
			pNA05:  5 km resolution 

	avn -		NOAA's aviation model
			avn2arl: 1 degree lat/lon input GRIB 

	coamps -	NAVY COAMPS model
			cmp2arl: version 2
			cmp3arl: version 3

	ecmwf -		ECMWF grib encoded files (hybrid or sigma)
			ecm2arl: requires GRIBEX libraries 

	eta -		NOAA's eta model
			eta04arl:  4 km
			eta12arl: 12 km 
			eta40arl: 40 AWIPS-212 grid

	gfs -		variation of grib2arl with more levels
			gfs2arl: consistent with current ARL global archive

	grib -		generic grib decoder for NOAA AVN & ECMWF
			grib2arl: general GRIB to ARL for pressure level GRIB 
			inventory: brief summary of GRIB records in file
			content: detailed summary of each GRIB record
			unpacker: sample program to unpack real data array

	mm5 -		converter from mm5/ver3 to ARL format
			mm5toarl:

	ncar -		decoder for NCAR/NCEP reanalysis files
			ncr2arl: 2.5 degree GRIB files

	netcdf -	NOAA's CDC netCDF reanalysis archive
			cdf2arl: 2.5 degree NetCDF data files

	nmm - 		NOAA's non-hydrostatic mesoscale model
			nmm2arl

	rams -		Regional Atmospheric Modeling System
			rams2arl.exe: converts RAMS history files

	rsm - 		NOAA's regional spectral model
			rsmp2arl: input data on pressure surfaces
			rsms2arl: input data on sigma surfaces

	sample - 	Sample program to create ARL packed files from
			user supplied data
			dat2arl

        wrfgrib2arl -   WRF output in GRIB format


Special Library directories:

	library	-	Master library directory that contains soft links
			to each of the libraries created in their respective
			directories. Each library must be compiled.

	libfcsubs -	"C" version of direct access I/O routines
			that are called by the various Fortran programs
			that read GRIB data files.

	libhysplit -	Selected subroutines from the main HYSPLIT library
			to pack and unpack ARL and GRIB formatted data 
			records.

	libcmapf -	Fortran programs to convert geographic positions
			between latitude-longitude coordinates and 
			conformal map coordinates.

