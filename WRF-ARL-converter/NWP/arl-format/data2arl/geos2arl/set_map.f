! Use emacs f90-mode:   -*- mode:f90  -*-
module set_map_module
  implicit none
  include 'netcdf.inc'
  private
  public :: set_map
contains
  subroutine set_map(diag,handle,nx,ny,grids,ierr)
    ! read lon,lat, construct grids info
    use setvar_module
    implicit none
    logical, intent(in) :: diag
    integer, intent(in) :: handle, nx, ny
    real, intent(out) :: grids(:)
    integer, intent(out) :: ierr

    real :: lon(nx), lat(ny)
    integer :: varid, n3d, ndim, dimlen(nf_max_var_dims)

    !  read in lon, lat
    call setvar(diag,handle,'lon',varid,n3d,ndim,dimlen,ierr)
    if (ierr .eq. 0) then
       if (ndim .ne. 1 .or. nx .ne. dimlen(1)) ierr=1
    end if
    if (ierr .ne. 0) then
       write (*,*) 'set_map: setvar error for lon'
       return
    end if
    ierr = nf_get_var_real(handle,varid,lon)
    if (ierr .ne. nf_noerr) then
       write (*,*) 'set_map: error getting lon values'
       return
    end if

    call setvar(diag,handle,'lat',varid,n3d,ndim,dimlen,ierr)
    if (ierr .eq. 0) then
       if (ndim .ne. 1 .or. ny .ne. dimlen(1)) ierr=1
    end if
    if (ierr .ne. 0) then
       write (*,*) 'set_map: setvar error for lat'
       return
    end if
    ierr = nf_get_var_real(handle,varid,lat)
    if (ierr .ne. nf_noerr) then
       write (*,*) 'set_map: error getting lat values'
       return
    end if

    ! Populate grids array for a lat-lon grid
!!$ In metset: reading grids(12) from extended header:
!!$   GRID(KG,KT)%POLE_LAT, GRID(KG,KT)%POLE_LON,  GRID(KG,KT)%REF_LAT,     &
!!$   GRID(KG,KT)%REF_LON,  GRID(KG,KT)%SIZE,      GRID(KG,KT)%ORIENT,      &
!!$   GRID(KG,KT)%TANG_LAT, GRID(KG,KT)%SYNC_XP,   GRID(KG,KT)%SYNC_YP,     &
!!$   GRID(KG,KT)%SYNC_LAT, GRID(KG,KT)%SYNC_LON,  GRID(KG,KT)%DUMMY,       &

!!$ In gblset: determine if this is a lat-lon grid
!!$  IF(GRID(KG,KT)%SIZE.EQ.0.0)THEN
!!$     GRID(KG,KT)%LATLON=.TRUE.

!!$ In gbl2ll: Grid system is simply defined as the number of grid points
!!$            from the corner point at 1,1 using an even lat-lon increment
!!$            x,y - grid position
!!$  CLAT=GRID(KG,KT)%SYNC_LAT+(Y-1.0)*GRID(KG,KT)%REF_LAT
!!$  CLON=GRID(KG,KT)%SYNC_LON+(X-1.0)*GRID(KG,KT)%REF_LON

    ! grids( 1: 2): GRID(KG,KT)%: POLE_LAT, POLE_LON: lat/lon of corner point 2 (nx,ny)
    grids(1) = lat(ny)
    grids(2) = lon(nx)
    ! grids( 3: 4): GRID(KG,KT)%: REF_LAT, REF_LON:   grid spacing in lat/lon
    grids(3) = lat(2)-lat(1)
    grids(4) = lon(2)-lon(1)
    ! grids( 5)   : GRID(KG,KT)%: SIZE: =0, flag for lat/lon grid
    grids(5) = 0.
    ! grids( 6: 9): GRID(KG,KT)%: ORIENT, TANG_LAT, SYNC_XP, SYNC_YP: not used
    grids(6:9) = 0.
    ! grids(10:11): GRID(KG,KT)%: SYNC_LAT, SYNC_LON: lat/lon of corner point 1 (1,1)
    grids(10) = lat(1)
    grids(11) = lon(1)
    ! grids(12)   : GRID(KG,KT)%: DUMMY: not used
    grids(12) = 0.

    return
  end subroutine set_map
end module set_map_module


