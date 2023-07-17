! Use emacs f90-mode:   -*- mode:f90  -*-
module interp_subs

  USE map_utils

  ! Module contains subroutines needed for interpolation of gridded data to lat/lon points
  ! map projection routines:
  !   cmap: all taken from hymodelc code
  !   wps: calling module_map_utils routines taken from WPS (and installed in stilt/src)

  ! $Id: interp_subs.f,v 1.11 2013/12/18 15:43:08 mellis Exp $
  ! Change Log:
  ! $Log: interp_subs.f,v $
  ! Revision 1.11  2013/12/18 15:43:08  mellis
  ! mods related to working with GEOS global grid
  !
  ! Revision 1.11  2013/12/18 15:43:08  mellis
  ! mods related to working with GEOS global grid
  !
  ! Revision 1.10  2013/12/17 18:25:21  trn
  ! Update global grid handling as in stilt/multi: gblset.f and hymodelc.f90
  !
  ! Revision 1.9  2013/07/31 19:45:14  trn
  ! First revisions for interp_receptors
  !
  ! Revision 1.8  2013/03/12 16:35:42  trn
  ! Additions to support extra-long label/header for large nlvl
  !
  ! Revision 1.7  2011/06/14 22:04:45  trn
  ! Mostly cosmetic changes; use zero buffer by default
  !
  ! Revision 1.6  2011/06/14 19:34:11  jhegarty
  ! Jennifer Hegarty: Added capability to include a buffer zone within domain
  !
  ! Revision 1.5  2008/06/13 16:52:58  trn
  ! Disabled module-scope implicit none statement
  !
  ! Revision 1.4  2008/06/13 13:12:21  trn
  ! Fine-tuned printout control
  !
  ! Revision 1.3  2008/01/17 18:40:43  trn
  ! Support use of WPS geolocation routines
  !
  ! Revision 1.2  2008/01/15 19:14:26  trn
  ! Fixed choice of mapping routine for WRF
  !
  ! Revision 1.1  2007/12/03 18:15:04  trn
  ! Added terrain interpolation program
  !

!  implicit none  !Interferes with implicit real, implicit integer statements in subroutines below for pgf90

  private

  CHARACTER(4), save               :: MODEL 
  integer, save ::          ICX,       MN,                                              &
         NX,       NY,        NZ,                                              &
         K_FLAG,   LENH
  real, save ::                                                   &
         POLE_LAT, POLE_LON,  REF_LAT,                                         &
         REF_LON,  SIZE,      ORIENT,                                          &
         TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
         SYNC_LAT, SYNC_LON,  DUMMY

  real, save :: gbase(15)

  TYPE(proj_info), save :: proj

  logical, save :: latlon, global, map_initialized=.FALSE.

  REAL, PARAMETER :: REARTH = 6371.2    ! radius of earth in km

  real, parameter :: RADPDG=pi/180., DGPRAD=180./pi, pi_2=pi/2.
  REAL, PARAMETER :: DEGPRD = 180./pi  ! deg per radian

  real, parameter :: vmiss=-999.99, vmissle=-999.

  public :: interp_init, interp_2d, interp_ll2xy

contains

  subroutine interp_init(nxin,nyin,header,maxlenh,ierr_interp,iprint,label)
    
    integer, intent(in) :: nxin,nyin,maxlenh,iprint
    character (len=maxlenh) :: header
    character(52) :: label
    character(2)  :: cgrid
    integer  :: knx, kny
    integer :: IYR,IMO,IDA,IHR,IFC
    character(4)  :: KVAR
    integer, intent(out) :: ierr_interp
    character(len=256) :: fmt

    integer :: ierr, k108, k50, recl_k108
    real :: clat1, clon1, clat2, clon2, dlat, dlon, DSX,DSY

    real :: dxm, reflon180, synclon180, true1, true2
    integer :: iproj


! decode extended portion of the header
    k108=108
    k50=50
    recl_k108 = k108+k50
    fmt='(A4,I3,I2,12F7.0,3I3,I2,I4)'
    if (LABEL(11:11) .eq. 'X') then
       k50=52
       recl_k108 = k108+k50
       write (*,'(a,i5)') &
            'NOTE: Detected extra-long label with label length k50=',k50
    endif
    if (header(105:105) .eq. 'X') then
       k108=110
       recl_k108 = k108+k50
       fmt='(A4,I3,I2,12F7.0,3I3,I2,1x,I5)'
       write (*,'(a,i5)') &
            'NOTE: Detected extra-long  header with standard header length k108=',k108
    end if

! From metset.f: read in cgrid
    if (k50 .eq. 50) then
       READ(label,'(5I2,2X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
    else
       READ(label,'(5I2,4X,a2,A4)')IYR,IMO,IDA,IHR,IFC,cgrid,KVAR
    end if
    
    READ(HEADER(1:k108),fmt,iostat=ierr)                &
         MODEL,    ICX,       MN,                                              &
         POLE_LAT, POLE_LON,  REF_LAT,                                         &
         REF_LON,  SIZE,      ORIENT,                                          &
         TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
         SYNC_LAT, SYNC_LON,  DUMMY,                                           &
         NX,       NY,        NZ,                                              &
         K_FLAG,   LENH

    if (ierr .ne. 0) then
       WRITE(*,'(a)') 'interp_init: ERROR: decoding header, cannot recover'
       WRITE(*,'(2a)') 'HEADER=',HEADER(1:k108)
       ierr_interp = 1
       return
    end if
 
    ! From metset.f:
    KNX=ICHAR(CGRID(1:1))
    KNY=ICHAR(CGRID(2:2))
    IF(KNX.GE.64.OR.KNY.GE.64)THEN
       NX=(KNX-64)*1000+NX
       NY=(KNY-64)*1000+NY
!!$     GRID(KG,KT)%NUMBER=KNX*10+KNY
!!$  ELSE
!!$     READ(CGRID,'(I2)')IGRID
!!$     GRID(KG,KT)%NUMBER=IGRID
    END IF
    
    if (nx .ne. nxin .or. ny .ne. nyin) write (*,*) 'interp_init: Warning: inconsistent nx, ny, nxin, nyin: ', &
         & nx, ny, nxin, nyin

    if (pole_lat .le. vmissle .and. pole_lon .le. vmissle .and. ref_lat .le. vmissle) then
       LATLON=.FALSE.
       GLOBAL=.FALSE.
       ! using WPS geolocation routines
       dxm = size*1000.
       true1 = tang_lat
       true2 = dummy
       iproj = nint(orient)
       synclon180 = sync_lon
       if (synclon180 > 180) synclon180=synclon180-360
       reflon180 = ref_lon
       if (reflon180 > 180) reflon180=reflon180-360
       call map_set(proj_code=iproj,proj=proj,lat1=sync_lat,lon1=synclon180,knowni=sync_xp,knownj=sync_yp, &
            & dx=dxm, stdlon=reflon180, truelat1=true1, truelat2=true2)
       gbase(:) = vmiss
       if (iprint .ge. 1) &
            & write(*,*) "interp_init: Using WPS geolocation routines with proj_code: ",iproj
    else

    !gblset:
       LATLON=.FALSE.
       GLOBAL=.FALSE.
       IF(SIZE .EQ. 0.0)THEN
          LATLON=.TRUE.
          CALL GBLDXY(1.0,1.0,DSX,DSY)
          SIZE=DSY
       ! find the corner points
          CLAT1=SYNC_LAT
          CLON1=SYNC_LON
          CLAT2=POLE_LAT
          CLON2=POLE_LON
       ! grid spacing
          DLAT=REF_LAT
          DLON=REF_LON
       ! determine if the grid is global
          IF((CLON2+DLON-CLON1 .GE. 360.0 .or. &
              CLON2+DLON-CLON1 .EQ.  0.0) .AND. &
              CLAT2-CLAT1 .GE. 180.0) GLOBAL=.TRUE.
       endif

       IF(MODEL.EQ.'RAMS')THEN
          if (iprint .ge. 1) WRITE(*,*) "interp_init: Using cmap code for RAMS with REF_LAT=",REF_LAT
          CALL SOBSTR(GBASE,                                 &
               &                  REF_LAT,REF_LON)
       ELSE
          if (iprint .ge. 1) WRITE(*,*) "interp_init: Using cmap code for not RAMS with TANG_LAT=",TANG_LAT
          IF(.NOT.LATLON)THEN
             CALL STLMBR(GBASE,                                 &
                  &                  TANG_LAT, REF_LON)
          END IF
       END IF
       IF(.NOT.LATLON)THEN
       !          use single point grid definition
          CALL STCM1P(GBASE,                                    &
               &               SYNC_XP,  SYNC_YP,               &
               &               SYNC_LAT, SYNC_LON,              &
               &               REF_LAT,  REF_LON,               &
               &               SIZE,     ORIENT)
       END IF
    endif
    if (iprint .ge. 1) then
       write(*,*) "interp_init: REF_LON: ",REF_LON
       write(*,*) "interp_init: SYNC_XP,SYNC_YP: ",SYNC_XP,SYNC_YP
       write(*,*) "interp_init: SYNC_LON,SYNC_LAT: ",SYNC_LON,SYNC_LAT
       write(*,*) "interp_init: SIZE,NX,NY: ",SIZE,NX,NY
    end if
    map_initialized = .TRUE.
    return
  end subroutine interp_init

  subroutine interp_2d(lon,lat,interr,nxin,nyin,header,maxlenh, &
       & outterr,ierr_interp,iprint,bufix,bufjy,interp_type,label)

    implicit none
    real, intent(in) :: lon, lat
    integer, intent(in) :: nxin,nyin,maxlenh,iprint,interp_type
    real, intent(in) :: interr(nxin,nyin)
    character (len=maxlenh) :: header
    character(52) :: label
    real, intent(out) :: outterr
    integer, intent(out) :: ierr_interp
    real, intent(in) :: bufix, bufjy

    real :: xp, yp, wt_i, wt_j
    integer :: i1, i2, j1, j2
    
    ierr_interp = 0

    if (.not. map_initialized) then

       call interp_init(nxin,nyin,header,maxlenh,ierr_interp,iprint,label)
       if (ierr_interp .ne. 0) return

    end if

    call interp_ll2xy (LAT,LON,XP,YP)

    IF(XP.LT.(1+bufix).OR.XP.GT.(NX-bufix) .OR.            &
         &      YP.LT.(1+bufjy).OR.YP.GT.(NY-bufjy)) THEN
       ierr_interp=2
       if(iprint .ge. 2) WRITE(*,*) 'point off grid: lat,lon= ',LAT,LON,' xp,yp= ',XP,  YP
       return
    END IF
    if (iprint .ge. 3) write (*,*) 'lat,lon= ',LAT,LON,' xp,yp= ',XP,  YP
    if (bufix .gt. 0 .and. iprint .ge. 3) write(*,*) '*** NX,1+bufix,nx-bufix,xp ***',NX,1+bufix,NX-bufix,XP
    if (bufjy .gt. 0 .and. iprint .ge. 3) write(*,*) '*** NY,1+bufjy,ny-bufjy,yp ***',NY,1+bufjy,NY-bufjy,YP
    
    if (interp_type .eq. 1) then
       ! bilinear interpolation
       i1 = int(xp)
       i2 = min(nx,i1+1)
       j1 = int(yp)
       j2 = min(ny,j1+1)
       
       wt_i=xp - real(i1)
       wt_j=yp - real(j1)
       outterr = (1.-wt_i)*(1.-wt_j)*interr(i1,j1) + &
            &     wt_i *(1.-wt_j)*interr(i2,j1) + &
            & (1.-wt_i)*    wt_j *interr(i1,j2) + &
            &     wt_i *    wt_j *interr(i2,j2)
       if (iprint .ge. 4) then
          write (*,*) 'i1,i2,wt_i= ',i1,i2,wt_i
          write (*,*) 'j1,j2,wt_j= ',j1,j2,wt_j
          write (*,*) 'll,lr,ul,ur= ',interr(i1,j1),interr(i2,j1),interr(i1,j2),interr(i2,j2)
          write (*,*) 'interp= ',outterr
       end if
    elseif (interp_type .eq. 0) then
       !nearest neighbor
       i1 = min(nx,nint(xp))
       j1 = min(ny,nint(yp))
       outterr = interr(i1,j1)
    else
       ierr_interp=3
       if(iprint .ge. 1) WRITE(*,*) 'interp_2d: invalid interp_type=',interp_type
    END IF
    return
  end subroutine interp_2d

  subroutine interp_ll2xy (LAT,LON,XP,YP)

    implicit none
    real, intent(in) :: lat, lon
    real, intent(out) :: xp, yp

    if (.not. map_initialized) stop 'interp_ll2xy: not map_initialized'

    if (all(gbase .le. vmissle)) then
       ! use WPS mapping routines
       call latlon_to_ij(proj, lat, lon, xp, yp)
    else
       ! use cmap mapping routines
       IF(LATLON)THEN
          CALL GBL2XY(LAT,LON,XP,YP)
          IF(GLOBAL)THEN
             IF(XP .GE. FLOAT(NX+1)) &
                  XP=XP-FLOAT(NX)
             IF(XP .LT. 1.0) XP=NX+XP
             IF(YP .GT. FLOAT(NY)) &
                  YP=2.0*NY-YP
             IF(YP .LT. 1.0) YP=2.0-YP
          end IF
       ELSE
          CALL CLL2XY(GBASE,LAT,LON,XP, YP)
       END IF
    endif

    return
  end subroutine interp_ll2xy

  SUBROUTINE GBL2XY(CLAT,CLON,X,Y)

    IMPLICIT NONE

    REAL,  INTENT(IN)  :: clat,clon                         ! latlon location
    REAL,  INTENT(OUT) :: x,y                               ! grid position

    REAL               :: tlat,tlon

! Grid system is simply defined as the number of grid points
! from the corner point at 1,1 using an even lat-lon increment
! for the x and y directions. Grid distances are computed
! where needed according to the latitude of the grid point

    TLAT=CLAT
    IF(TLAT.GT. 90.0)TLAT= 180.0-TLAT
    IF(TLAT.LT.-90.0)TLAT=-180.0-TLAT
    Y=1.0+(TLAT-SYNC_LAT)/REF_LAT

    TLON=CLON
    IF(TLON.LT.0.0)  TLON=360.0+TLON
    IF(TLON.GT.360.0)TLON=TLON-360.0
    TLON=TLON-SYNC_LON
    IF(TLON.LT.0.0)TLON=TLON+360.0
    X=1.0+TLON/REF_LON

    return

  end SUBROUTINE GBL2XY

  subroutine gbldxy(x,y,gsx,gsy)
    implicit none
    real, intent(in) :: x,y
    real, intent(out) :: gsx, gsy
    real :: clat, clon

    CALL GBL2LL(X,Y,CLAT,CLON)

    ! latitude grid spacing
    GSY = (2.0*PI*REARTH) / (360.0/REF_LAT)

    ! longitude grid spacing
    ! GSX = COS(CLAT/DEGPRD)*GSY*REF_LON/REF_LAT
    GSX = COS(CLAT/DEGPRD)*GSY*REF_LON/REF_LAT

    return
  end subroutine gbldxy

  subroutine gbl2ll(x,y,clat,clon)
    implicit none
    REAL,    INTENT(IN)  :: x,y         ! grid position
    REAL,    INTENT(OUT) :: clat,clon   ! latlon location

    CLAT=SYNC_LAT+(Y-1.0)*REF_LAT
    IF(CLAT.GT. 90.0)CLAT= 180.0-CLAT
    IF(CLAT.LT.-90.0)CLAT=-180.0-CLAT

    CLON=SYNC_LON+(X-1.0)*REF_LON
    CLON=MOD(CLON,360.0)
    IF(CLON.GT.180.0)CLON=CLON-360.0

    return
  end subroutine gbl2ll

  subroutine sobstr(stcprm, p_lat, p_lon)
    !*
    !*  Set Map Parameters for an OBlique STeReographic Projection
    !*  Inputs: p_lat,P_lon - latitude and longitude of the
    !*          projection pole (tangent point)
    !*  Outputs: stcprm - map parameters
    !*/
    !
    ! Original RCS Id: sobstr.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    dimension stcprm(15)
    TEMP=90.
    ! CHG(09/15/03) pass on pole lat and lon (temp2 and temp4) opposing ref_lat and lon
    ! might need adjustment for northern hemisphere
    TEMP2=-(180.0+p_lat)
    TEMP4=p_lon+180.0
    !      call mpstrt(stcprm, TEMP, p_lat,p_lon, TEMP2,p_lon)  Original: NAN's
    ! CHG(09/15/03) pass on pole lat and lon (temp2 and temp4) opposing ref_lat and lon
    !      WRITE(*,*)'mpstrargs:',TEMP,TEMP2,TEMP4, p_lat,p_lon
    call mpstrt(stcprm, TEMP, TEMP2,TEMP4, p_lat,p_lon)
    return
  END subroutine sobstr

  subroutine mpstrt(stcprm, conang, p_lat, p_long, r_lat, r_long)
    !
    ! Original RCS Id: mpstrt.f90,v 1.2 2006/06/19 15:30:42 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    !*  REARTH=6356.766  from U.S. Standard Atmosphere, 1976
    !*  REARTH=6367.47   for spherical earths in GRIB grids
    !*  REARTH=6371.2    original assumption for CMAPF routines.
    !*                   source lost, probably old NWPC grids.

    ! CHG(09/15/03) adapt to saulos version
    !      PARAMETER (REARTH=6367.47)
    PARAMETER (REARTH=6367.00)
    dimension stcprm(15)
    !typedef struct {
    ! 1:     gamma,
    ! 2-10:  rotate(3,3),
    ! 11-12: x0,y0,
    ! 13-14: crotate,srotate,
    !  15:   grdszq
    !

    !*
    !*  General Purpose routine to Set Map Parameters.  Called by the
    !*  special purpose routines for oblique stereographic, oblique and
    !*  transverse Mercator, oblique Lambert Conformal, etc. Projections
    !*  Inputs: p_lat,p_long - Latitude and Longitude of the Pole
    !*            Point for the Projection
    !*          r_lat,r_long - Latitude and Longitude of a
    !*            reference point - 180degrees around the Pole Point
    !*            from the cut
    !*          conang - angle between the Projection Pole axis and
    !*            the generators (sides) of the cone on which the Earth
    !*            is projected. + or - 90 degrees indicates a plane,
    !*            hence Stereographic; 0 degrees indicates a cylinder
    !*            and hence Mercator.
    !*  Outputs: stcprm - map parameters
    !*/
    dimension temp(3)
!    REAL, EXTERNAL :: x_prod !now inside this module, not external

    ifind(l,k) = -2 + 3*l + k !statement function

    call ll_geo(p_lat,p_long,temp)
    do k=1,3
       stcprm(ifind(3,k)) = temp(k)
    enddo
    !  for (k=0;k<3;k++) stcprm->rotate[2][k] = temp.v[k];
    call ll_geo(r_lat,r_long,temp)
    do k=1,3
       stcprm(ifind(1,k)) = temp(k)
    enddo
    !  for (k=0;k<3;k++) stcprm->rotate[0][k] = temp.v[k];
    tnorm = x_prod(stcprm(ifind(3,1):ifind(3,1)+2),stcprm(ifind(1,1):ifind(1,1)+2), &
         &                 stcprm(ifind(2,1):ifind(2,1)+2))
    !  norm =
    !  x_product (stcprm->rotate[2],stcprm->rotate[0],stcprm->rotate[1]);
    do k=1,3
       stcprm(ifind(2,k)) = stcprm(ifind(2,k)) /tnorm
    enddo
    !  for (k=0;k<3;k++) stcprm->rotate[1][k] /= norm
    tnorm = x_prod(stcprm(ifind(2,1):ifind(2,1)+2),stcprm(ifind(3,1):ifind(3,1)+2), &
         &                 stcprm(ifind(1,1):ifind(1,1)+2))
    !  x_product (stcprm->rotate[1],stcprm->rotate[2],stcprm->rotate[0]);
    !  stcprm->x0 = stcprm->y0 = stcprm->srotate = 0;
    stcprm(11)=0.
    stcprm(12)=0.
    stcprm(14)=0.
    stcprm(13)=1.
    stcprm(15) = REARTH
    stcprm(1) = sin(RADPDG * conang )
    ! stcprm->crotate = 1.;
    ! stcprm->gridszeq = REARTH;
    ! stcprm->gamma = sin(RADPDEG * cone_ang);
    !*geographic triple : i = equator @ std merid,
    !       j = equator @ 90E,
    !       k = North Pole */
    !*map base triple : i' = M-prime meridian @ M-equator,
    !     j' = M-East at location i',
    !     k' = M-Pole */
    !
    return
  end subroutine  MPSTRT

  subroutine ll_geo(xlat, xlong, vector)
    !  Given a latitude xlat and longitude xlong, returns the unit vector
    !  directed to the given point in the geo (geographic) coordinate system.
    !
    !  The geo system is a 3-D Cartesian coordinate system, with origin
    !  at the Center of the Earth, with the x_3 axis pointing to the North
    !  Pole, the x_1 axis pointing to the Equator at the Greenwich Meridian,
    !  and the x_2 axis by the right hand rule pointing to the Equator
    !  at 90 E.
    !
    ! Original RCS Id: ll_geo.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    dimension vector(3)
    vector(3) = sin(RADPDG * xlat)
    clat = cos(RADPDG * xlat)
    vector(1) = clat * cos(RADPDG * xlong)
    vector(2) = clat * sin(RADPDG * xlong)
    return
  END subroutine ll_geo

  subroutine stlmbr(stcprm, reflat, reflon)
    !*
    !*  Set Map Parameters for a North Polar LaMBeRt Conic Conformal
    !*    Projection
    !*  Inputs: reflat - tangent latitude of the tangent latitude of
    !*                   cone.
    !*          reflon - midrange longitude (180 degrees from cut)
    !*  Outputs: stcprm - map parameters
    !*/
    !
    ! Original RCS Id: stlmbr.f90,v 1.1 2005/09/27 13:39:34 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    TEMP=90.
    TEMP2=0.
    CALL mpstrt(stcprm,reflat,TEMP,TEMP2,reflat,reflon)
    return
  END subroutine stlmbr

  subroutine stcmap(stcprm, tnglat, reflon)
    !*
    !*  Set Map Parameters for a North Polar LaMBeRt Conic Conformal
    !*    Projection
    !*  included for compatibliity with previous version
    !*  Inputs: tnglat - tangent latitude of the tangent latitude of
    !*                   cone.
    !*          reflon - midrange longitude (180 degrees from cut)
    !*  Outputs: stcprm - map parameters
    !*/
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    TEMP=90.
    TEMP2=0.
    CALL mpstrt(stcprm,tnglat,TEMP,TEMP2,tnglat,reflon)
    return
  END subroutine stcmap

  subroutine stcm1p(stcprm, x1, y1, xlat1, xlong1,                  &
       & xlatrf, xlonrf, gridsz, orient)
    !
    ! Original RCS Id: stcm1p.f90,v 1.2 2006/02/28 16:20:16 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    real enx,eny,enz,norm,x1a,y1a
!    REAL, EXTERNAL :: cgszll !now inside this module, not external

    c_or = cos(RADPDG * orient)
    s_or = - sin(RADPDG * orient)
    !  stcprm->x0 = stcprm->y0 = stcprm->srotate = 0;
    stcprm(11) = 0.
    stcprm(12) = 0.
    stcprm(13) = 1.
    stcprm(14) = 0.
    stcprm(15) = 1.
    !  stcprm->crotate = stcprm -> gridszeq = 1.0;
    call cpolll(stcprm, xlatrf, xlonrf, enx, eny, enz)
    norm = sqrt(enx*enx + eny*eny)
    if (norm .eq. 0.) then
       call cgrnll(stcprm,xlatrf,xlonrf,enx,eny,enz)
       norm = sqrt (enx* enx + eny*eny)
    endif
    enx = enx/norm
    eny = eny/norm
    !  stcprm->gridszeq *= gridsz / cgszll(stcprm, xlatrf,xlonrf);
    stcprm(15) = stcprm(15) * gridsz / cgszll(stcprm,xlatrf,xlonrf)

    !  stcprm -> crotate = eny * c_or - enx * s_or;
    stcprm(13) = eny * c_or - enx * s_or
    !  stcprm -> srotate = -eny * s_or - enx * c_or;
    stcprm(14) = -eny * s_or - enx * c_or
    call cll2xy(stcprm, xlat1,xlong1, x1a,y1a)

    !  stcprm->x0 += x1 - x1a;
    !  stcprm->y0 += y1 - y1a;
    stcprm(11) = stcprm(11) + x1 - x1a
    stcprm(12) = stcprm(12) + y1 - y1a
    return
  END subroutine stcm1p

  subroutine cpolll(stcprm, alat,along, enx,eny,enz)
    ! returns a vector aligned with the direction toward the North Pole.  I.e.
    ! parallel to the vector from the Earth's center to the North Pole.
    !
    ! Original RCS Id: cpolll.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    real map(3),pole(3),geog(3)
    call ll_geo(alat,along, geog)
    call basg2m(stcprm, geog, map)
    do k=1,3
       pole(k) = stcprm(3*k + 1)
    enddo
    call proj_3d(stcprm, map, pole, enx,eny,enz)
    return
  END subroutine cpolll

  subroutine cgrnll(stcprm, alat,along, enx,eny,enz)
    ! returns a vector aligned in the direction toward the Greenwich Meridian
    ! at the equator.  I.e. parallel to the vector from Earth's center to
    ! that point.
    !
    ! Original RCS Id: cgrnll.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    real map(3),pole(3),geog(3)
    call ll_geo(alat,along, geog)
    call basg2m(stcprm, geog, map)
    do k=1,3
       pole(k) = stcprm(3*k - 1)
    enddo
    call proj_3d(stcprm, map, pole, enx,eny,enz)
    return
  END subroutine cgrnll

  subroutine cll2xy(stcprm, alat, along, x, y )
    !
    ! Original RCS Id: cll2xy.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    real geog(3),temp(3)
    call ll_geo(alat,along, geog)
    call basg2m(stcprm, geog, temp)
    call map_xy(stcprm, temp, x,y)
    return
  END subroutine cll2xy

  subroutine basg2m(stcprm, xin, xout)
    !  receives the vector xin, giving the components of a point in the geo
    !  system; returns xout, the components of the same point in the map system.
    !
    !  Conversions to and from the geo (geographic) and the map
    !  (Map-oriented)  coordinate systems.
    !
    !  The geo system is a 3-D Cartesian coordinate system, with origin
    !  at the Center of the Earth, with the x_3 axis pointing to the North
    !  Pole, the x_1 axis pointing to the Equator at the Greenwich Meridian,
    !  and the x_2 axis by the right hand rule pointing to the Equator
    !  at 90 E.
    !
    !  In the map system, the axis of the map projection passes through the
    !  center of the Earth, and through a point on the surface that acts as the
    !  "Pole" of the map (Which will coincide with the North or South Pole,
    !  unless the projection is Oblique or Transverse).  The x_3 axis of the
    !  map coordinate system is aligned with this map pole.  In Lambert and
    !  Mercator projections, a "Cut" extends from the map pole along a great
    !  circle to its antipode; the unrolled map is discontinuous there.  The
    !  x_1 axis of the map coordinate system is diametrically opposite this
    !  cut, and 90 degrees from the pole; the x_2 axis is selected to complete
    !  a right hand rule.
    !
    !  The coefficients of the map coordinate system relative to the geo system
    !  are given in elements 2 through 10 of the stcprm array.
    !
    ! Original RCS Id: basg2m.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    dimension stcprm(15),xin(3),xout(3)
    ifind(l,k) = -2 + 3*l + k
    do k=1,3
       xout(k) = 0.
    enddo
    do l=1,3
       do k=1,3
          xout(l) = xout(l) + xin(k) * stcprm(ifind(l,k))
       enddo
    enddo
    return
  END subroutine basg2m

  subroutine proj_3d(stcprm, point, vect, enx,eny,enz)
    !/*
    ! *  At a given point, resolves the components of a vector vect in the
    ! *  local coordinate system of the map projection.  It is assumed that
    ! *  vect and point is given in the 3-dimensional geocentric _map_
    ! *  coordinate system, rather than the North centered _geo_ system.
    ! *  returns these components as enx, eny, and enz.  Replaces vect with
    ! *  its projection on the tangent plane at point.
    ! */
    !
    ! Original RCS Id: proj_3d.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real stcprm(15)
    real point(3),vect(3)
    dot_pr = 0.
    do k=1,3
       dot_pr = dot_pr + point(k)*vect(k)
    enddo
    !/*
    ! *  dot_prod is the local vertical component of vect.  Note, point is
    ! *  assumed to be a unit vector.
    ! */
    enz = dot_pr
    do k=1,3
       vect(k) = vect(k) - dot_pr * point(k)
    enddo
    !/*  vect has now been projected to a plane tangent to point.
    ! */
    fact = 1. + point(3)
    xi = vect(1) - vect(3) * point(1) / fact
    eta = vect(2) - vect(3) * point(2) / fact
    !/*
    ! *  xi, eta represent components of vector on the projection plane,
    ! *  i.e. the 2-dimensional map system.
    ! */
    fact = stcprm(1) -1.
    if (fact .lt. 0.) then
       !/*
       ! *  For Lambert Conformal and Mercator projections (gamma < 1.0) ,
       ! *  a rotation of the vector components is needed.
       ! */
       if ( abs(point(3)) .lt. 1.) then
          theta = fact * atan2(point(2),point(1))
          cgthta = cos(theta)
          sgthta = sin(theta)
          fact = xi * cgthta - eta * sgthta
          eta = eta * cgthta + xi * sgthta
          xi = fact
       endif
    endif
    !/*
    ! *  Now rotate xi, eta to the final map projection direction.
    ! */
    enx = eta * stcprm(13) - xi * stcprm(14)
    eny = - xi * stcprm(13) - eta * stcprm(14)
    return
  END subroutine proj_3d

  subroutine map_xy(stcprm, x_map, x, y)
    !
    ! Original RCS Id: map_xy.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    dimension stcprm(15),x_map(3)
    call map_xe(stcprm, x_map, xi, eta, 'c')
    call xe_xy(stcprm, xi, eta, x, y)
    return
  END subroutine map_xy

  subroutine map_xe(stcprm, x_map, xi, eta,flag)
    !
    ! Original RCS Id: map_xe.f90,v 1.1 2005/09/27 13:39:33 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    character * 1 flag
    dimension stcprm(15),x_map(3)
    if (abs(x_map(3)) .ge. 1.) then
       !* Projection pole or its antipodes
       xi = 0.
       if (stcprm(1) * x_map(3) .gt. 0.) then
          eta = 1./stcprm(1)
       else
          if (x_map(3) .gt. 0) then
             eta = 1.0e10
          else
             eta = -1.0e10
          endif
       endif
       return
    else
       if (abs(stcprm(1)) .eq. 1.) then
          !* Stereographic Case, away from pole
          fact = 1. / (stcprm(1) + x_map(3) )
          xi = x_map(2)*fact
          eta = 1./stcprm(1) - x_map(1) * fact
       else
          ymerc = .5 * log( ( 1. + x_map(3) ) / (1. - x_map(3) ) )
          !* This Projection has a cut.  Select theta according to the rule
          !* If cutflag = 'e', -PI/2 <= theta < 3PI/2, if cutflag = 'E',
          !* 0 <= theta < 2PI, if cutflag = 'w', -3Pi/2 <= theta < PI/2,
          !* cutflag = 'W', -2PI <= theta < 0., else -PI<=theta<PI.
          if (flag .eq. 'E') then
             theta = atan2(- x_map(2), - x_map(1) ) + pi
          else if (flag .eq. 'W') then
             theta = atan2(- x_map(2), - x_map(1) ) - pi
          else if (flag .eq. 'e') then
             theta = atan2(- x_map(1), x_map(2) )  + pi_2
          else if (flag .eq. 'w') then
             theta = atan2(x_map(1), - x_map(2) )  - pi_2
          else
             theta = atan2(x_map(2), x_map(1) )
          endif
          rhog = xpabva( stcprm(1), - ymerc)
          !* rhog = ( exp( - gamma * ymerc ) - 1. ) / gamma */
          xi = (1. + stcprm(1) * rhog) * snabva(stcprm(1), theta)
          eta = stcprm(1) * (1. + stcprm(1) * rhog) *                 &
               &                        csabva(stcprm(1), theta) - rhog
       end if
    end if
    return
  END subroutine map_xe

  subroutine xe_xy(stcprm, xi, eta, x, y)
    !
    ! Original RCS Id: xe_xy.f90,v 1.1 2005/09/27 13:39:34 trn Exp $
    !
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    PARAMETER (REARTH=6367.47)
    dimension stcprm(15)
    x = stcprm(11) + REARTH / stcprm(15) *                          &
         &  (stcprm(13) * xi + stcprm(14) * eta)
    y = stcprm(12) + REARTH / stcprm(15) *                          &
         &  (stcprm(13) * eta - stcprm(14) * xi)
    return
  END subroutine xe_xy

  real function csabva(a, b)
!
! $Id: interp_subs.f,v 1.11 2013/12/18 15:43:08 mellis Exp $
!
    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

!* returns (1. - cos(a * b) ) / a / a, or its limit as a -> 0
    term = snabva(.5 * a, b)
    csabva = .5 * term * term
    return
  end function csabva
      
  real function snabva(a, b)
!* returns sin(a * b) / a, or its limit as a -> 0.

    IMPLICIT REAL (A-H,O-Z)
    IMPLICIT integer (I-N)

    real a,b,c,csq
    c = a * b
    csq = c * c
    if (csq .gt. .001) then
       snabva = sin( c )  / a
    else
       snabva = b * ( 1. - csq / 6. *                                &
            &                 ( 1. - csq / 20. *                               &
            &                 ( 1. - csq / 42. )))
    endif
    return
  end function snabva
      
!  A series of functions which may have problems in certain ranges
!  of one of their parameters.  These functions, if written normally,
!  may suffer round-off problems.  They are used in particular for
!  Lambert Conformal projections whose parameters make them almost
!  the same as a Mercator Projection.

  real function xpabva(a, b)
!* returns (exp (a * b) - 1. ) / a, or its limit as a -> 0.
    real a,b,c,csq
    c = a * b
    csq = .25 * c * c
    if ( csq  .gt. .001) then
       xpabva = (exp( c ) - 1.) / a
    else
       xpabva = exp(.5 * c) * b * (1. + csq / 6. *                  &
            &                               (1. + csq / 20. *                  &
            &                               (1. + csq / 42. )))
    endif
    return
  end function xpabva

  real FUNCTION x_prod(A,B,C)
!  Returns as C the vector cross product of A and B:
!  C = A x B.  Also returns as function value the absolute value of C.
!
! $Id: interp_subs.f,v 1.11 2013/12/18 15:43:08 mellis Exp $
!
    IMPLICIT NONE

    REAL, INTENT(IN)  :: a(3),b(3)
    REAL, INTENT(OUT) :: C(3)

    INTEGER, PARAMETER :: indx(4) = (/2,3,1,2/)
    INTEGER            :: k
    REAL    :: D

    D = 0.
    do k=1,3
       C(k) = A(indx(k))*B(indx(k+1)) - A(indx(k+1))*B(indx(k))
       D = D + C(k) * C(k)
    enddo
    x_prod = SQRT(D)

  END FUNCTION x_prod

  real function cgszll(stcprm, alat, along)
!
! $Id: interp_subs.f,v 1.11 2013/12/18 15:43:08 mellis Exp $
!
    IMPLICIT NONE
    real stcprm(15), alat, along

    real map(3)
    real geog(3)
    real ymerc

    call ll_geo(alat, along, geog)
    call basg2m(stcprm, geog, map)
    if (map(3) .ge. 1.) then
       if (stcprm(1) .ge. 1.) then
          cgszll = 2. * stcprm(15)
       else
          cgszll = 0.
       endif
    else if (map(3) .le. -1.) then
       if (stcprm(1) .le. -1.) then
          cgszll = 2. * stcprm(15)
       else
          cgszll = 0.
       endif
    else if (abs(stcprm(1)) .ge. 1.) then
       cgszll = stcprm(15) * (1. + sign(map(3),stcprm(1)))
    else
       ymerc = -.5 * log( (1. - map(3))/(1. + map(3)) )
       cgszll = stcprm(15) * exp( - (1.-stcprm(1)) * ymerc) *         &
            &   (1. + map(3))
    endif
    return
  end function cgszll
  
end module interp_subs
