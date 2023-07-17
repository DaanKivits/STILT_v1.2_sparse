./data2arl/avn  			Last Revised: 17 Nov 2003 
----------------------------------------------------------------------
This directory contains the old grib decoder for NOAA's Aviation model
grib files with data on pressure surfaces. The input GRIB files are
assumed to be one degree resolution.  Output options are a default
100x100 grid point conformal subgrid, a hemispheric polar stero-
graphic projection, or a global latitude-longitude grid. This decoder
has been replaced by the more general (and complicated) program in 
the grib directory, which will handle any resolution input global
model GRIB output file from NOAA or ECMWF data sources.
