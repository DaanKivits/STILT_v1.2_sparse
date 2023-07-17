GRIB1 file processing with data on a regular latitude-longitude grid

-------------------------------------------------------------------------

 Usage: grib2arl [-options]
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
  -r[rain fall accumulation time hours: {6}]
  -t[the number of time periods to process: {744}]
  -z[zero initialization of output file 0:no {1}:yes]
  
 Note-1: for some Linux systems {export MALLOC_CHECK_=0}
