! Use emacs f90-mode:   -*- mode:f90  -*-
module setmap_module

! $Id: setmap.f,v 1.4 2010/09/24 19:07:32 trn Exp $

contains
SUBROUTINE setmap (diag,ncid,nxp,nyp,grids,ierr,geoloc_type)

  IMPLICIT NONE

  include 'netcdf.inc'

  LOGICAL,       INTENT(IN)  :: diag
  INTEGER,       INTENT(IN)  :: ncid
  INTEGER,       INTENT(IN)  :: nxp,nyp
  character (len=4) :: geoloc_type
  integer,       INTENT(out) :: ierr
  real,          INTENT(out) :: grids(12)

  CHARACTER(80) :: label
  INTEGER       :: status,iproj
  REAL          :: clat,clon,orient,plat,tlat,tlon,dxkm,dykm
  real          :: polar_lat, polar_lon, tang_lat

  real, parameter :: vmiss=-999.99

  polar_lat=90.
  polar_lon=0.

  ierr=0
  label='MAP_PROJ'
  status = NF_GET_ATT_INT(ncid, nf_global, label, iproj)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 1
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',iproj
  END IF

  label='CEN_LAT'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, clat)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 2
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',clat
  END IF

  label='CEN_LON'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, clon)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 3
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',clon
  END IF

  label='STAND_LON'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, orient)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 4
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',orient
  END IF

  label='DX'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, dxkm)
  dxkm=dxkm/1000.0
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 5
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',dxkm
  END IF

  label='DY'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, dykm)
  dykm=dykm/1000.0
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 6
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',dykm
  END IF

  label='TRUELAT1'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, tlat)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 7
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',tlat
  END IF

  label='TRUELAT2'
  status = NF_GET_ATT_REAL(ncid, nf_global, label, plat)
  IF(status.NE.nf_noerr) THEN
     WRITE(*,*) trim(label),':',NF_STRERROR(status)
     ierr = 8
     return
  ELSE
     IF(diag)WRITE(*,*)trim(label),' = ',plat
  END IF

!--------------------------------------------------------
! HYSPLIT packed data configuration

! grid orientation
  GRIDS(6)=0.0     
! delta=x grid size in km
  GRIDS(5)=0.5*(dxkm+dykm)
! synch point in x,y coordintes
  GRIDS(8)=0.5*(NXP+1)
  GRIDS(9)=0.5*(NYP+1)
! synch point in lat/lon coordinates
  GRIDS(10)=CLAT 
  GRIDS(11)=CLON 

! variable reserved for future use
  GRIDS(12)=0.0

! defines a polar sterographic projection
  IF(iproj.EQ.2)THEN

!    set the pole position and reference lat/lon
     label='POLE_LAT'
     status = NF_GET_ATT_REAL(ncid, nf_global, label, polar_lat)
     IF(status.NE.nf_noerr) THEN
        WRITE(*,*) trim(label),':',NF_STRERROR(status)
        ierr = 8
        return
     ELSE
        IF(diag)WRITE(*,*)trim(label),' = ',polar_lat
     END IF
     label='POLE_LON'
     status = NF_GET_ATT_REAL(ncid, nf_global, label, polar_lon)
     IF(status.NE.nf_noerr) THEN
        WRITE(*,*) trim(label),':',NF_STRERROR(status)
        ierr = 8
        return
     ELSE
        IF(diag)WRITE(*,*)trim(label),' = ',polar_lon
     END IF
     GRIDS(1)=polar_lat
!    pole longtitude (+180 from cut)
     GRIDS(2)=polar_lon

!    reference lat/lon (at which grid size specified)
     GRIDS(3)=TLAT
!    reference longitude and grid alignment
     GRIDS(4)=ORIENT

!    tangent latitude
     GRIDS(7)=TLAT 

! defines a mercator projection
  ELSEIF(iproj.EQ.3)THEN

!    pole lat/lon axis through pole
     GRIDS(1)=0.0
     GRIDS(2)=ORIENT

!    reference lat
     GRIDS(3)=TLAT
!    reference lon
     GRIDS(4)=ORIENT

!    tangent latitude
     GRIDS(7)=TLAT

! defines a lambert conformal projection
  ELSEIF(iproj.EQ.1)THEN

! Added for consistency with GRIB specification
     if (tlat .lt. 0) polar_lat = -90.
!    pole lat/lon axis through pole
     GRIDS(1)=polar_lat  !not: TLAT
     GRIDS(2)=polar_lon  !not: ORIENT

!    reference lat
     tang_lat=eqvlat(tlat,plat)
     GRIDS(3)=Tang_LAT
!    reference lon
     GRIDS(4)=ORIENT

!    tangent latitude
     GRIDS(7)=Tang_LAT

  ELSE
     WRITE(*,*)'Undefined projection: ',iproj
     ierr = 9
  END IF
  
  if (trim(geoloc_type) .eq. 'wps' .or. trim(geoloc_type) .eq. 'WPS') then
! Redefine parts of the header for the WPS geolocation routines
     grids(1) = vmiss
     grids(2) = vmiss
     grids(3) = vmiss
     grids(6) = real(iproj)
     if (iproj == 1) then
        grids(7) = tlat
        grids(12) = plat
     endif
  else
     if (trim(geoloc_type) .ne. 'cmap' .and. trim(geoloc_type) .ne. 'CMAP') &
          & write (*,*) 'Unsupported geoloc_type=',trim(geoloc_type),' using default=cmap'
  end if

  return

END SUBROUTINE setmap 
      real function eqvlat(lat1,lat2)
!*  Written 12/21/94 by Dr. Albion Taylor
!*
!*    This function is provided to assist in finding the tangent latitude
!*    equivalent to the 2-reference latitude specification in the legend
!*    of most lambert conformal maps.  If the map specifies "scale
!*    1:xxxxx true at 40N and 60N", then eqvlat(40.,60.) will return the
!*    equivalent tangent latitude.
!*  INPUTS:
!*    lat1, lat2:  The two latitudes specified in the map legend
!*  RETURNS:
!*    the equivalent tangent latitude
!*  EXAMPLE:  stcmap(& strcmp, eqvlat(40.,60.), 90.)
!*/
!/*   Changes Made May 9, 2003 to accomodate the following special
!*   situations:
!*   1. if lat1 == lat2, returned value will be lat1 (reduced to between -90.
!*      and 90.).
!*   2. If either lat1 or lat2 is 90. (or -90.) then 90. (or -90.) will be
!*      returned.  This reflects the fact that, for y fixed (-90. < y < 90.),
!*      as x -> 90. (or -90.), eqvlat(x,y) ->90. (or -90.)  This limiting
!*      case of tangent latitude 90. is a polar stereographic projection,
!*      for which the scale at 90. is a maximum, and therefore greater than the
!*      other latitude y.  Thus, eqvlat(90.,60.) returns 90., although the
!*      scale at 90. will be greater than at 60. For eqvlat(90.,-90.), the
!*      limit process is ambiguous; for the sake of symmetry, such a case
!*      will return 0.
!*/

        implicit none

!beg:     include 'cmapf.fi'
!*  REARTH=6356.766  from U.S. Standard Atmosphere, 1976
!*  REARTH=6367.47   for spherical earths in GRIB grids
!*  REARTH=6371.2    original assumption for CMAPF routines.
!*                   source lost, probably old NWPC grids.

      real pi, pi_2, dgprad, radpdg
      PARAMETER (pi=3.14159265358979,pi_2=pi/2.)
      parameter (dgprad=180./pi,radpdg=pi/180.)

      real rearth
      PARAMETER (REARTH=6367.47)
      integer k_gama, k_reflon, k_x0, k_y0
      parameter (k_gama=1, k_reflon=k_gama+1)
      parameter (k_x0=k_reflon+1, k_y0=k_x0+1)
      integer k_crot, k_srot, k_gdszeq, k_n1dgr, k_s1dgr, k_maparam
      parameter (k_crot=k_y0+1, k_srot=k_crot+1, k_gdszeq=k_srot+1)
      parameter (k_n1dgr=k_gdszeq+1,k_s1dgr=k_n1dgr+1)
      parameter (k_maparam=k_s1dgr)
!end:     include 'cmapf.fi'

      real lat1,lat2

      real fsm, slat1, slat2, temp, tau, al1, al2
      parameter (fsm = 1.e-3)
        slat1 = sin(radpdg * lat1)
        slat2 = sin(radpdg * lat2)
! reorder slat1, slat2
        if (  slat1 .lt. slat2) then
          temp = slat1
          slat1 = slat2
          slat2 = temp
        endif
!/*  Take care of special cases first */
        if (slat1 .eq. slat2) then
          eqvlat = asin(slat1) * DGPRAD
          return
        endif
        if (slat1 .eq. -slat2 ) then
          eqvlat = 0.
          return
        endif
        if (slat1 .ge. 1.) then
          eqvlat = 90.
          return
        endif
        if (slat2 .le. -1.) then
          eqvlat = -90.
          return
        endif
!/* Compute al1 = log((1. - slat1)/(1. - slat2))/(slat1 - slat2) */
        tau = (slat1 - slat2)/(2. - slat1 - slat2)
        if ( tau .gt. FSM )then
          al1 = log((1. - slat1)/(1. - slat2)) / (slat1 - slat2)
        else
          tau = tau * tau
          al1 = -2./(2. - slat1 - slat2) * &
     &                  (1.    + tau * &
     &                  (1./3. + tau * &
     &                  (1./5. + tau * &
     &                  (1./7. + tau))))
        endif
!/* Compute al2 = log((1. + slat1)/(1. + slat2))/(slat1 - slat2) */
        tau = (slat1 - slat2)/(2. + slat1 + slat2)
        if ( tau .gt. FSM ) then
          al2 = log((1. + slat1)/(1. + slat2)) / (slat1 - slat2)
        else
          tau = tau * tau
          al2 =  2./(2. + slat1 + slat2) * &
     &                  (1.    + tau * &
     &                  (1./3. + tau * &
     &                  (1./5. + tau * &
     &                  (1./7. + tau))))
        endif
        eqvlat = asin ((al1 + al2) / (al1 - al2)) * DGPRAD
        return
      end function eqvlat

end module setmap_module
