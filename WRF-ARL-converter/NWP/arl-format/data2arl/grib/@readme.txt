./data2arl/grib				Last Revised: 18 Mar 2004  
----------------------------------------------------------------------

SYNOPSIS

This directory contains a generic grib decoder for ECMWF and NOAA
formatted grib files output from the global model in latitude-
longitude grid coordinates. The decoder processes one or multiple
time periods per execution. Data may be on pressure, sigma, or
the hybrid vertical coordinate.  Only the minimum number of fields
required to run hysplit are processed.

Convert Latitude-Longitude GRIB data file to ARL format with the data
organized by time so that a consecutive group of records contains all
the variables at the same time.  The ARL format consists of an index
record followed by data records.  One record per variable per level,
then followed by records for the next time period.  All records are of
fixed length and packed one byte per variable.  Packing information is
coded in the header portion of each record.

The program repacks the input data into ARL format in two different
modes. If a latitude-longitude output grid is defined (-g3 or -g4) then
the latlon input grid is repacked and written directly to ARL format.
A latlon input grid may be resized as a subgrid by setting the number of
output grid points (-n{nxp:nyp}).

The latlon input data may also be interpolated to a conformal map
projection (-g0, -g1, -g2) by specifying the output grid size (km),
grid center latlon, and the number of grid points in each direction.
The projection is automatically defined depending upon the latitude
of the grid center, where Polar Sterographic is used for latitudes
above 60, Mercator for latitudes less than 30, and Lambert Conformal
for intermediate latitudes. The extract conformal map always defaults
to a 100 km resolution projection of 100x100 grid points centered at
a lat/lon point set on the command line.

----------------------------------------------------------------------
DATA AVAILABILITY

NOAA global model grib files may be obtained from ftpprd.ncep.noaa.gov
These files are packaged with all the fields together.

ECMWF global grib files (ERA-40) are available from:
http://data.ecmwf.int/data, where the grib files are constructed by '
the user from a selection of fields.  It is recommended that at a
minimum the following fields be extracted:

Pressure Levels: geopotential, temperature, u-velocity, v-velocity,
		vertical velocity, relative humidity

Surface: 2 m temperature, 10 m U & V wind components, total precipitation

Invariant: geopotential

----------------------------------------------------------------------

USAGE: grib2arl [-options]

  -i[primary grib data: file name {required}]
  -s[supplemental grib data: file name {optional}]
  -c[constant grib data: file name {optional}]
  -x[subgrid extract center longitude {-80.0}]
  -y[subgrid extract center latitude {60.0}]
  -g[output projection  0 :conformal extract
                        1 :fixed northern hemisphere polar
                        2 :fixed southern hemisphere polar
                       {3}:lat-lon global grid (as input)
                        4 :lat-lon extract grid
  -n[number of (x:y) extract grid points {100}]
  -k[number of output levels including sfc {16}]
  -p[surface defined by {1}:pressure or 0:terrain height]
  -q[analyze grib file {0} or use saved configuration: 1]
  -z[zero initialization of output file 0:no {1}:yes]

-i: Specifies the primary input grib file that contains the 3D
    meteorological fields. Data may be for one time period or 
    multiple time periods in the same input grib file. File
    may also contain surface fields, but everything should be
    time ordered.

-s: Specifies the supplemental grib file that contains additional
    fields that may not be in the primary input file such as
    precipitation, 10 m winds, and 2 m temperatures. If multiple
    time periods are processed, then this file should contain
    exactly the same number of time periods as the primary file.

-c: Constant time invariant fields, such as geopotential, which
    is converted to terrain height by the decoder.  To properly
    use the output file, it must contain either terrain  height
    or surface pressure.  Frequently these fields are not contained
    in the 2d or 3d grib files.

-x: Center longitude of the output grid if an extraction grid is 
    specified as option -g0 or -g4. The default extraction grid 
    resolution is 100 km for -g0.

-y: Center latitude of the output grid if an extraction grid is 
    specified as option -g0 or -g4.

-g: Specifies the output grid. Options 0,1,2 result in the input 
    being interpolated to a conformal map projection.  Options 1
    and 2 are predefined northern and southern hemispheric polar
    stereographic projections.  Option 3 is the no-interpolation
    case, where the output data are written to the same latitude
    longitude global grid as the input data. Option 4 defines an
    extract of the latitude-longitude grid. In the case of -g4,
    the input grid is required to be global. Otherwise use -g3.

-n: Specifies the number of grid points (both x and y) for the
    extraction grid option -g0 or -g4. Different dimensions for 
    x and y can be specified -n140:70 defines a 140x by 70y grid.

-k: Specifies the number of levels (including the surface) counted
    from the ground in the output grid.  This number may have to
    be adjusted after several trials if the input data contain
    fewer levels than the default value (16).  Grib records do
    not contain any information regarding how many levels are
    available in the file.
  
-p: Specifies the use of surface pressure (1) or terrain height (0)
    as one of the output surface fields. Model calculations with
    the output data require some information about the surface
    elevations.  The default is to use surface pressure.  One of
    fields should be available in the input grid data. This option
    is only valid when both fields are available in the input data. 

-q: Will read the CFG_GRIB file if available. This causes the 
    program to skip the data set analysis step and use the input
    grib file configuration created from the last time the program
    was run.  This can save considerable time when processing
    identical multiple files.

-z: Causes the output field to be initialized prior to converting    
    processing the input grib data.  The default is initialization.
    This should be turned off, if for instance one wants to use the
    forecast fields from the previous time period, a new execution
    of grib2arl, will then write the new fields into the old file,
    only replacing those fields that are available.  An example
    might be precipitation, that would only be available from the
    previous forecast and not the model initialization file.

----------------------------------------------------------------------

GRIB2ARL OUTPUT FILES:

DATA.ARL	- (binary) the final processed data in ARL format
CFG_ARL	   	- (ascii) the configuration of the packed data
CFG_GRIB        - (binary) describes the input grib file configuration
MESSAGE		- (ascii) processing diagnostic and error messages

Note that multiple time period output files can be constructed from
multiple executions of grib2arl with single time period input files
by using the unix "cat" command:

cat DATA{time2}.ARL >> DATA{time1}.ARL

----------------------------------------------------------------------

Several other utility routines are provided that can be used as examples
to decode grib data records.  The following routines DO NOT use any 
w3lib routines and are entirely self contained.

	inventory - dumps short listing of records in grib file
	content - complete listing of each section in a grib record
	unpacker - decodes grib real data array 

