!----------------------------------------------------------------
! University of Houston
! Air Quality Modeling and Monitoring Center
! 218 Old Science Building
! Houston, TX 77204-5048
! USA
!
! email: sbkim@math.uh.edu (Dr. Seung-Bum Kim)
!        dwbyun@math.uh.edu (Dr. Daewon W. Byun)
!----------------------------------------------------------------
!
!-----RAMS fields
real, dimension(mnx) :: xrm, xrt
real, dimension(mny) :: yrm, yrt  
real, dimension(mnx,mny) :: xrutm, yrutm

!-----ARL and intermediate ARRAY
! 3D
real, dimension(mnz,mnx,mny) ::        &
     u,v,w,t,p,rv,dn0,tke,pi0,         &
     pp,theta              
! 2D
real, dimension(mnx,mny) ::            &
     umof,vmof,                        &
     dlwf,dswf,ustr,tstr,qstr,         &
     shgt,glat,glon,tpp1,              &
     rh2m,prss,t02m,                   &
     uwnd,vwnd,sphu,temp,pres,tken,wwnd 

! COMMON BLOCKS
common /ramsgrid/ xrm,yrm,xrt,yrt,xrutm,yrutm
common /arlfld/ u,v,w,t,p,rv,dn0,tke,pi0,pp,       &
                theta,umof,vmof,                   &
		dlwf,dswf,ustr,tstr,qstr,shgt,     &
		glat,glon,tpp1,rh2m,prss,t02m,     &
                uwnd,vwnd,sphu,temp,pres,tken,wwnd
