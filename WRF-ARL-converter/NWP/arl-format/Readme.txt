This directory contains the codes and related files for converting and
using ARL format files.

CVS tags:

 Branches:
   ARL: import branch for data2arl codes obtained from NOAA ARL tar files

 Other tags:
   import-start: NOAA ARL release from tar files obtained 2005/06/27 from:
     http://www.arl.noaa.gov/ready/hysp_data2arl.html.  The separate
     tar files were combined using the following procedure:
       emperor[122]> ls downloads/
       Readme.txt	  date		  hysp_data2arl.html  profile.zip
       afwa.tar.Z.tar	  display.zip	  mm5.tar.Z.tar       rams.tar.Z.tar
       arl2grad.zip	  eta.tar.Z.tar   ncar.tar.Z.tar      rsm.tar.Z.tar
       avn.tar.Z.tar	  gfs.tar.Z.tar   netcdf.tar.Z.tar    sample.tar.Z.tar
       coamps.tar.Z.tar  grib.tar.Z.tar  nmm.tar.Z.tar       showgrid.zip
       emperor[123]> foreach tarfile ( downloads/*tar.Z.tar )
       foreach> gunzip -c $tarfile | tar xvf - 
       foreach> end
     NOTE: the import ignored the library subdirectory and its
           contained symbolic links:
	     emperor[167]> ll $DAPD/p1229-carbon/arl-format/data2arl/library/
	     total 5
	     ./
	     ../
	     libcmapf.a -> ../libcmapf/libcmapf.a
	     libfcsubs.a -> ../libfcsubs/libxlffcsubs.a
	     libhysplit.a -> ../libhysplit/libhysplit.a

   Release-1-0: Version as uploaded to DEAS computer on 2006/01/16
   Release-1-3: Version as provided as part of WRF-STILT-Release-1-3, Sep 2007
   Release-2-1: Version as provided as part of WRF-STILT-Release-2-1, Sep 2007
   Release-2-2: Version as provided as part of WRF-STILT-Release-2-2, Oct 2007
   Release-2-3: Version as provided as part of WRF-STILT-Release-2-3, Feb 2008
   pleiades_v32_01, pleiades_v32_02, pleiades_v32_03: Version ported to Pleiades, Apr 2009
   Release-3-2-1: Version as provided as part of WRF-STILT-Release-3-2-1, Oct 2010

Directories and their main executables:

   interp_terrain.x: program to interpolate data from arl file to lat/lon
   arl-chk.csh and chk_data.x: script and program to check arl-format files
   data2arl: reformatting codes
     wrfgrib2arl: WRF GRIB -> ARL reformatter (deprecated, untested for WRFV3)
     wrfnc2arl: WRF netcdf -> ARL reformatter

To import a new NOAA ARL release:

   untar the distribution file, then cd into data2arl

   import the S/W into the appropriate repository, using:
   - cvs import -m"release comment" NWP/arl-format/data2arl ARL release-tag
  
To build the libraries and executables:
  edit and execute install-arl.csh
  NOTE: invoke without arguments for help message
        Make sure to issue an approriate "setenv STILTSRC ..." command
        if interp_terrain.x is required
      
To run the wrfnc2arl executable:
   1. Make a local copy of data2arl/wrfnc2arl/var_sample and edit as needed
   2. Run "data2arl/wrfnc2arl/wrfnc2arl -P var_sample fname"
   3. Rename DATA_*.WRF and CFG_*.wrf to names unique to fname
   4. Repeat steps 2 and 3 for all output times

      Sample c-shell commands, with variables set as follows:
         wrfdir - directory containing the wrfout... files to be processed
         arldir - directory containing data2arl directory tree of
                  codes
	 foreach gfile ( $wrfdir/wrfout* )
	  $arldir/data2arl/wrfnc2arl/wrfnc2arl \
	     -P var_sample \
	     $gfile >&! $gfile:t.out
	  foreach dfile ( DATA_*.WRF CFG_*.WRF )
	   /bin/mv $dfile $gfile:t.$dfile
	  end
	 end
   
Other files:

   Readme.txt - this file
   chk_data.f - standalone F90 program for reading ARL file and
                printing out sample values
   chk_data.in - input file for chk_data.f
   split_chk_data.csh - reformats output from chk_data.f for subsequent
                        analysis by:
   read.split.chkdata.q - Splus code for ingest and plotting of
                        split_chk_data.csh output
   read.split.chkdata.r - R code for ingest and plotting of
                        split_chk_data.csh output
