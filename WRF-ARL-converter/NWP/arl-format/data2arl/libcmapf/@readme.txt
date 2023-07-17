The Cmapf routines form a library of FORTRAN-language subroutines and functions
which process geographic and dynamic data between grids and geographical
coordinates on a family of conformal map projections.  Data characterizing the
map projections and the grids on them is embedded in maparam structures. Note:
in all cases, North latitudes will be given as positive degrees, South latitudes
negative degrees, East longitudes positive degrees and West longitudes negative
degrees.  Author and questions: albion.taylor@noaa.gov


Initializing (configuring) the maparam structure:
-------------------------------------------------

Stage I - Selecting a Map Projection.

To initialize a map parameter (maparam) structure with a Conformal Polar
Projection;

dimension stcprm(9)

call stlmbr(stcprm, tnglat, reflon)

will fill the structure "stcprm" with default data specifying a Lambert
Conformal Projection on a cone tangent to the globe at latitude tnglat degrees.
If a cut is needed to lay the cone out flat on the plane, it will be made
opposite to (180 degrees away from) longitude reflon.  If tnglat = 90 or -90,
the map becomes a North or South Polar Stereographic Projection (and the value
of reflon is immaterial).  If tnglat=0., the map becomes a Mercator Projection.

For non-zero values of tnglat, between (but not including) -90 and +90, the
structure defines a Lambert Conformal projection tangent to the Earth at
"tnglat".  If the desired tangent latitude is not known, but the map legend
specifies two reference latitudes, the tangent latitude may be found by

tnglat = eqvlat(reflat1,reflat2)

Caution: this value is not the same as the mean latitude 0.5*(reflat1+reflat2).
Attempts to use the mean latitude will yield misleading results.

Stage II - Laying out a grid on a selected map projection.

This stage configures the maparam structure with the information to locate the
points of the grid on the selected map projection.  Two alternative means are
provided: Either One-Point specification or Two-Point specification.

a) Two-Point Specification:

call stcm2p(stcprm, x1,y1,xlat1,xlon1, x2,y2,xlat2,xlon2)

causes stcprm to be configured so that the grid point with grid coordinates
x1,y1 will be located on the map at latitude xlat1, longitude xlon1, while
x2,y2 is located on the map at xlat2, xlon2.  The two points must be different.
This specification locates all other grid points, assumed arranged on the
projection with equal spacing in the x- and y- directions.

The two points will often be opposite corners of the intended grid, but many
times a different selection will be preferred.  For example, in a Mercator
projection, it will often be better to select two points of the same latitude or
the same longitude, to ensure exact N-S orientation of the grid.  Alternatively,
several subgrids of an overall grid may be guaranteed compatibility by using the
same two points of the overall grid, even if neither of the points lie in the
subgrid at all.

b) One-Point Specification:

call stcm1p(stcprm, x1,y1,xlat1,xlon1, scalat,scalon, gsize,orient)

causes stcprm to be configured so that the grid point with grid coordinates
x1,y1 will be located on the map at latitude xlat1, longitude xlon1.  The
remaining configuration information is provided at the scaling point whose
latitude and longitude are scalat, scalon.  At that point, and incidentally
every on the latitude circle scalat, the size of a grid cell is gsize in both
the x- and y- directions.  Finally, at the scaling point, the y-lines of the
grid make an angle of "orient" degrees clockwise (East) of the scalon meridian.

The x1,y1 coordinates can be thought of as "pinning" that grid point to the map
at xlat1,xlon1, the gsize value can be thought of as "zooming" the grid in and
out until the desired size is attained, and the "orient" specification can be
thought of as rotating the grid about the pinned point until the grid has the
wanted orientation.

Frequently, a grid will be published as having an "lov" value of a longitude
intended to be vertical; for such grids, scalon=lov and orient=0.

Using the configured maparams:
-----------------------------

Coordinate transforms:
----------------------

dimension stcprm(9)

call cll2xy(stcprm, xlat,xlon, x,y)

will return grid coordinates x,y according to the grid configured as stcprm,
for the point whose latitude and longitude are xlat,xlon.  Longitude values are
adjusted to fit between reflon-180. and reflon+180., so that nearby points on the
globe remain nearby points on the map, unless they straddle the "cut"  at
longitude reflon+180., and find themselves on opposite edges of an unrolled
cylinder.

call cxy2ll(stcprm, x,y, xlat,xlon)

will return latitude and longitude parameters according to the grid configured
as stcprm, for the point whose grid coordinates are x,y.  Longitude values will
be returned in the range from -180. to +180.
If the point x,y corresponds to the North or South Pole, where longitude is
undefined, the returned longitude willl be 180. for the North Pole, 0. for the
South Pole, for compatibility with the WMO standard for reporting winds there.
(In an earlier version of CMAPF, the value of reflon was returned here.)

Wind Vector transforms:
-----------------------

dimension stcprm(9)

call cc2gll(stcprm, xlat,xlon, ue,vn, ug,vg)
and
call cg2cll(stcprm, xlat,xlon, ug,vg, ue,vn)

will convert the East- and North- components ue,vn (Cartographic components) of
the wind to and from components ug,vg of the wind in the x- and y- directions of
the grid (Grid components).  The first transfers from cartographic components to
grid components, the second does the reverse. In both cases, the wind is that
observed at the specified cartographic coordinates xlat,xlon.
At the poles (xlat=90. or xlat= -90.), the direction of "North", which
determines the meaning of ue and vn, is based on the nominal value of longitude.
The direction of North will agree with that perceived by an observer stationed a
short distance away on the xlon meridian.  Thus, if xlat,xlon = 90.,-70.,
"North" will be asssigned to the direction of the 110E meridian.

call cc2gxy(stcprm, x,y, ue,vn, ug,vg)
and
call cg2cxy(stcprm, x,y, ug,vg, ue,vn)

will convert the East- and North- components ue,vn (Cartographic components) of
the wind to and from components ug,vg of the wind in the x- and y- directions of
the grid (Grid components).  The first transfers from cartographic components to
grid components, the second does the reverse. In both cases, the wind is that
observed at the specified grid coordinates x,y.
If x,y refers to the North or to the South pole, the orientation of "North" is
defined according to World Meteorological Organization (WMO) code 878.  The
direction of "North" is that of the Greenwich meridian, which is the equivalent
of using xlat=90.,xlon=180. for the North Pole, and xlat=-90.,xlon=0. for the
South Pole.


call cw2gll(stcprm, xlat,xlon, ue,vn, ug,vg)
call cg2wll(stcprm, xlat,xlon, ug,vg, ue,vn)
call cw2gxy(stcprm, x,y, ue,vn, ug,vg)
and
call cg2wxy(stcprm, x,y, ug,vg, ue,vn)

perform the same services as the above transforms, except that within 1 degree
(111km or 60nm) or either pole, ue and vn are interpreted according to WMO code
878 above.  In an earlier version of CMAPF, this was the definition of the
c2g and g2c routines, above.

Gridsize Evaluation:
--------------------

dimension stcprm(9)

gsize = cgszll(stcprm, xlat,xlon)
and
gsize = cgszxy(stcprm, x,y)

return the size, in kilometers, of a grid cell at the indicated point.

Curvature Vector:
-----------------
dimension stcprm(9)

call ccrvxy (stcprm, x,y, gx,gy)

call ccrvll (stcprm, xlat,xlon, gx,gy)

return, in gx and gy, the logarithmic gradient of the scale.  This vector
represents the curvature on the Earth of a straight line segment on the map.
Units are in radians per kilometer on the Earth's surface.  See the paper for
details on usage in equations.

North Polar Vector evaluations:
-------------------------------
dimension stcprm(9)

call cpolxy (stcprm, x,y, enx,eny,enz)

call cpolll (stcprm, xlat,xlon, enx,eny,enz)


return, in enx, eny, and enz the components of the direction of the Earth's
axis.  components enx and eny provide the local direction of North in grid
coordinates, while enz is proportional to the local Coriolis parameter.
