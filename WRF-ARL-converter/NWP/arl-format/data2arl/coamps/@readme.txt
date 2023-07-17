./data2arl/coamps                                         Revised: 25 Apr 2002
------------------------------------------------------------------------------
Contains program to read multi-file Navy COAMPS model ouput
files and generate a single ARL packed meteorological data file.

cmp2arl - converts version 2 format coamps data files
cmp3arl - converts version 3 format coamps data files

There are several differences between version 2 and 3.  In particular the 
file name convention has changed and 3d variables are now written by level
in each file rather than writing a single 3d variable level to one file.
The header file (datahd_sfc....) has been expanded from 500 variables. Note
no documentation exists describing the new format and the grid description
information has changed.  Outputs appear to be interpolated to either standard
pressure levels or height levels.  Not all variables are available for each
vertical output type.  Vertical motion fields are not available although they
could be computed by integrating the divergence field (only available on
pressure level output).  Data fields are no longer written IEEE unformatted,
(4 byte header and trailer) but require direct access with the record length 
equaling 4 times the number of grid points. 

Code Notes:

1.  We modified your cmp2arl.f90, not your cmp3arl.f90, and include it
here as cmp3narl.f90.

2.  All changes are noted with either a !JJS or !PFC.

3.  We ran this on an pentium-linux machine, and therefore needed to 
add the endian-swap subroutine called SWAP32.

Problems:

1.  There is currently an error in the COAMPS 3 header data - the
longitude of the 1,1 nest grid point is incorrect (set to the XLON1
variable in SUB HEADER)in the current version of COAMPS.  We have
notified the COAMPS people at NRL-Monterey for a fix, but I am unsure
when it will be implemented, as it has no effect on any COAMPS data 
and is only included in the output header.  To workaround this, we
hard-coded these data points for our runs (see lines 222-231).

2.  We used a lambert conformal projection with the projection
intersection at 30 and 60 deg.  Therefore, in SUB WINDEX, we changed
GRIDS(3) for MTYPE.EQ.2 to equal REFLT2 (30 deg), a variable which had
to be added to the subroutine call list.  GRIDS(4) was left set equal
to ALIGN.  The other projections are untested.

3.  We processed hours 0-11 in each data set and therefore changed the
DO limit on line 114.

-------------------------------------
Peter F. Caffrey, Ph.D.
Naval Research Laboratory, Code 7228
4555 Overlook Ave SW
Washington, DC 20375

(202) 767-8474
(202) 404-8011 (fax)

Peter.Caffrey@nrl.navy.mil
