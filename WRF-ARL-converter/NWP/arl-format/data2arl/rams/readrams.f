!---------------------------------------------------------------
! University of Houston
! Air Quality Modeling and Monitoring Center
! 218 Old Science Building
! Houston, TX 77204-5048
! USA
!
! email: sbkim@math.uh.edu (Dr. Seung-Bum Kim)
!        dwbyun@math.uh.edu (Dr. Daewon W. Byun)
!----------------------------------------------------------------

subroutine readrams(nx,ny,nz,nzm1,lfirst,ng,flnm)
!----------------------------------------------------------------
! Name:     readrams
! Purpose:  read RAMS output files (analysis files) 
! Revised:  06 Aug 2002  Created for RAMS2ARL (S.-B. Kim)
! Note:     this version is F90 for use with RAMS v4.3
!
!---------------------------------------------------------------
!
include 'param.inc' 
include 'fields.inc' 
logical lfirst
character(len=80) :: flnm
character(len=32) :: cdname,cdunits
integer, parameter :: nind=300

real, allocatable :: au(:),av(:),aw(:),atempk(:),apress(:),      &
     avapor(:),adn0(:),atke(:),api0(:),app(:),        &
     atheta(:),                           &  !/3d
     auw(:),avw(:),arlong(:),arshort(:),           &
     austar(:),atstar(:),arstar(:),             &
     atopo(:),alat(:),alon(:),aprecip(:)              !/2d
real, allocatable :: ascratch(:)
integer, parameter :: mnxyz=mnx*mny*mnz
!
!--------------------------------------------------------------------
! Allocate necessary variables
!--------------------------------------------------------------------
allocate ( au(mnxyz) )
allocate ( av(mnxyz) )
allocate ( aw(mnxyz) )
allocate ( atempk(mnxyz) )
allocate ( apress(mnxyz) )
allocate ( avapor(mnxyz) )
allocate ( adn0(mnxyz) )
allocate ( atke(mnxyz) )
allocate ( api0(mnxyz) )
allocate ( app(mnxyz) )
allocate ( atheta(mnxyz) )
allocate ( auw(mnxyz) )
allocate ( avw(mnxyz) )
allocate ( arlong(mnxyz) )
allocate ( arshort(mnxyz) )
allocate ( austar(mnxyz) )
allocate ( atstar(mnxyz) )
allocate ( arstar(mnxyz) )
allocate ( atopo(mnxyz) )
allocate ( alat(mnxyz) )
allocate ( alon(mnxyz) )
allocate ( aprecip(mnxyz) )
allocate ( ascratch(mnxyz) )

!--------------------------------------------------------------------
!################## 3D FIELDS
!--->U: x-direction wind components [m/s]
call RAMS_varlib('u',nx,ny,nz,ng,AU,ascratch,&
      flnm,cdname,cdunits,nind) 
!--->V: y-direction wind components [m/s]
call RAMS_varlib('v',nx,ny,nz,ng,AV,ascratch,&
      flnm,cdname,cdunits,nind) 
!--->W: z-direction wind components [m/s]
call RAMS_varlib('w',nx,ny,nz,ng,AW,ascratch,&
      flnm,cdname,cdunits,nind) 
 
!--->TEMPK: temperature [K]
call RAMS_varlib('tempk',nx,ny,nz,ng,ATEMPK,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->PRESS: pressure [mb]
call RAMS_varlib('press',nx,ny,nz,ng,APRESS,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->VAPOR: water vapor mixing ratio [g/kg]
call RAMS_varlib('vapor',nx,ny,nz,ng,AVAPOR,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->DN0: reference state density [kg/m3]
call RAMS_varlib('dn0',nx,ny,nz,ng,ADN0,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->TKE: turbulent kinetic energy [m**2/s**2]
call RAMS_varlib('tke',nx,ny,nz,ng,ATKE,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->PI0: reference Exner function [J/(kgK)]
call RAMS_varlib('pi0',nx,ny,nz,ng,API0,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->PERT_PRESSURE: perturbation pressure [mb]
call RAMS_varlib('pert_pressure',nx,ny,nz,ng,APP,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->THETA: potential temperature [K]
call RAMS_varlib('theta',nx,ny,nz,ng,ATHETA,ascratch,&
      flnm,cdname,cdunits,nind) 

!################## 2D FIELDS

!--->UW [m2/s2]
call RAMS_varlib('uw',nx,ny,1,ng,AUW,ascratch,&
      flnm,cdname,cdunits,nind) 
 
!--->VW [m2/s2]
call RAMS_varlib('vw',nx,ny,1,ng,AVW,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->RSHORT [W/m2]
call RAMS_varlib('rshort',nx,ny,1,ng,ARSHORT,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->RLONG [W/m2]
call RAMS_varlib('rlong',nx,ny,1,ng,ARLONG,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->PRECIP (mm) --> converts to [m] & designate one of TPP1, TPP3....
call RAMS_varlib('precip',nx,ny,1,ng,APRECIP,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->USTAR_PS [m/s]
call RAMS_varlib('ustar_ps',nx,ny,1,ng,AUSTAR,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->TSTAR_PS [K]
call RAMS_varlib('tstar_ps',nx,ny,1,ng,ATSTAR,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->RSTAR_PS [kg/kg]
call RAMS_varlib('rstar_ps',nx,ny,1,ng,ARSTAR,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->TOPO [m]      
call RAMS_varlib('topo',nx,ny,1,ng,ATOPO,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->LAT [deg]
call RAMS_varlib('lat',nx,ny,1,ng,ALAT,ascratch,&
      flnm,cdname,cdunits,nind) 

!--->LON [deg]
call RAMS_varlib('lon',nx,ny,1,ng,ALON,ascratch,&
      flnm,cdname,cdunits,nind) 

!-----------------------------------------------------------------
!-----Load individual arrays from scratch vector a()
!
! --- 3D variables
do k = 1,nz 
  do j = 1,ny 
  do i = 1,nx 
    n3d = i + (j-1)*nx + (k-1)*nx*ny
    u(k,i,j) = au(n3d) 
    v(k,i,j) = av(n3d) 
    w(k,i,j) = aw(n3d) 
    t(k,i,j) = atempk(n3d) 
    p(k,i,j) = apress(n3d) 
    rv(k,i,j) = avapor(n3d)/1.e3  ! [g/kg] to [kg/kg] 
    dn0(k,i,j) = adn0(n3d) 
    tke(k,i,j) = atke(n3d) 
    pi0(k,i,j) = api0(n3d) 
    pp(k,i,j) = app(n3d) 
    theta(k,i,j) = atheta(n3d) 
  enddo 
  enddo 
enddo
! --- 2D variables
do k = 1,1 
  do j = 1,ny 
  do i = 1,nx 
    n3d = i + (j-1)*nx + (k-1)*nx*ny
    umof(i,j)=auw(n3d)
    vmof(i,j)=avw(n3d)
    dlwf(i,j)=arlong(n3d)
    dswf(i,j)=arshort(n3d)
    tpp1(i,j)=aprecip(n3d)/100. ! [mm] to [m]
    ustr(i,j)=austar(n3d)
    tstr(i,j)=atstar(n3d)
    qstr(i,j)=arstar(n3d)
    shgt(i,j)=atopo(n3d)
    glat(i,j)=alat(n3d)
    glon(i,j)=alon(n3d)
  enddo 
  enddo 
enddo
!
! Deallocate necessary variables
!
deallocate (au,av,aw,atempk,apress)
deallocate (avapor,adn0,atke,api0,app)  
deallocate (atheta)
deallocate (auw,avw,arlong,arshort)
deallocate (austar,atstar,arstar)
deallocate (atopo,alat,alon,aprecip) 
deallocate (ascratch)
!
return
end
