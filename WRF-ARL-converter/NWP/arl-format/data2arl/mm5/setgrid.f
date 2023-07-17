SUBROUTINE setgrid (ix, jx, kx, nicg, njcg, proj, gnum, dist, clat, clon, &
                    tru1, tru2, pole, cgi1, cgj1, cgi2, cgj2,             &
                    ka, n0, n1, vchar0, vchar1, arlsig)

!-------------------------------------------------------------------------------
! Name:     Setgrid      
! Purpose:  Defines mapping of latitude-longitude to grid units and initializes
!           ARL packing subroutines 
! Notes:    Uses standard cmapff and packer libraries 
! Revised:  14 Aug 2000  Original version.  (RRD)
!           09 Mar 2004  Carry through domain number (RRD)
!-------------------------------------------------------------------------------

      USE file

      IMPLICIT NONE

      INTEGER,       INTENT(IN)    :: ix
      INTEGER,       INTENT(IN)    :: jx
      INTEGER,       INTENT(IN)    :: kx
      INTEGER,       INTENT(IN)    :: nicg  
      INTEGER,       INTENT(IN)    :: njcg 
      INTEGER,       INTENT(IN)    :: proj
      INTEGER,       INTENT(IN)    :: gnum
      REAL,          INTENT(IN)    :: dist   
      REAL,          INTENT(IN)    :: clat 
      REAL,          INTENT(IN)    :: clon
      REAL,          INTENT(IN)    :: tru1 
      REAL,          INTENT(IN)    :: tru2 
      REAL,          INTENT(IN)    :: pole 
      REAL,          INTENT(IN)    :: cgi1 
      REAL,          INTENT(IN)    :: cgj1 
      REAL,          INTENT(IN)    :: cgi2 
      REAL,          INTENT(IN)    :: cgj2 
      INTEGER,       INTENT(IN)    :: ka,n0,n1
      CHARACTER*4,   INTENT(IN)    :: vchar0  ( : )
      CHARACTER*4,   INTENT(IN)    :: vchar1  ( : )
      REAL,          INTENT(IN)    :: arlsig  ( : )

      REAL                         :: parmap  ( 9 )
      REAL                         :: grids   (12 )
      CHARACTER*1                  :: num
      CHARACTER*80                 :: fname
      INTEGER                      :: k,n
      REAL,          EXTERNAL      :: eqvlat, cgszll
      REAL                         :: clat1,clon1,clat2,clon2,clatx,clonx
  

!-------------------------------------------------------------------------------
! Initialize latitude-longitude conversion for the coarse grid
! Note that in ARL system i and j conventions opposite those in MM5
!-------------------------------------------------------------------------------
  
  CALL stlmbr(parmap,eqvlat(tru1,tru2),clon)
  CALL stcm1p(parmap,0.5*(1.0+njcg),0.5*(1.0+nicg),clat,clon, &
              eqvlat(tru1,tru2),clon,dist,0.0)

!-------------------------------------------------------------------------------
! Compute latitude-longitude of dot domain corners and then reinitialize routine
! for the domain under consideration.
!-------------------------------------------------------------------------------

  CALL cxy2ll(parmap,cgj1,cgi1,clat1,clon1)
  CALL cxy2ll(parmap,cgj2,cgi2,clat2,clon2)
  CALL stcm2p(parmap,1.,1.,clat1,clon1,FLOAT(jx),FLOAT(ix),clat2,clon2)

! compute grid center position for diagnostic evaluation
  CALL cxy2ll(parmap,0.5*(1.0+jx),0.5*(1.0+ix),clatx,clonx)
  WRITE(*,'(a,2f10.2)')'Center grid lat-lon: ',clatx,clonx

!-------------------------------------------------------------------------------
! Create array that will be written to initialize ARL data packing
!-------------------------------------------------------------------------------

  IF( proj == 1 )THEN   ! lambert conformal

     grids(1) = pole               ! pole position
     grids(2) = clon               ! 180 from cut longitude
     grids(3) = eqvlat(tru1,tru2)  ! lat of grid spacing 
     grids(4) = clon               ! lon of grid spacing
     grids(5) = cgszll(parmap,grids(3),clon)  ! grid size
     grids(7) = eqvlat(tru1,tru2)  ! cone angle

  ELSEIF( proj == 2 )THEN ! polar sterographic

     grids(1) = pole               ! pole position
     grids(2) = 0.0                ! 180 from cut longitude
     grids(3) = eqvlat(tru1,tru2)  ! lat of grid spacing 
     grids(4) = clon               ! lon of grid spacing
     grids(5) = cgszll(parmap,grids(3),clon)  ! grid size
     grids(7) = 90.0               ! cone angle

  ELSEIF( proj == 3 )THEN ! mercator projection 

     grids(1) = pole               ! pole position
     grids(2) = 0.0                ! 180 from cut longitude
     grids(3) = eqvlat(tru1,tru2)  ! lat of grid spacing 
     grids(4) = clon               ! lon of grid spacing
     grids(5) = cgszll(parmap,grids(3),clon)  ! grid size
     grids(7) = 0.0                ! cone angle

  END IF 

  grids(6)  = 0.0                  ! orientation
  grids(8)  = 1.0                  ! lat/lon of grid 1,1
  grids(9)  = 1.0    
  grids(10) = clat1 
  grids(11) = clon1 
  grids(12) = 0.0                  ! unused

!-------------------------------------------------------------------------------
! Write ARL data packing configuration file                      
!-------------------------------------------------------------------------------

  write(num,'(I1)') (gnum/10)
  fname = 'CFG'//NUM//'_MM5'
  OPEN(iutarl,FILE=fname)      
  WRITE(iutarl,'(20X,a3,a1)')'MM5',NUM       ! abbreviated data set name
  WRITE(iutarl,'(20X,i4)') gnum,1            ! id and coordinate system
  WRITE(iutarl,'(20X,f10.2)')grids           ! orientation array
  WRITE(iutarl,'(20X,i4)')jx,ix,ka+1         ! domain and levels (+1 for sfc)
 
! surface level information
     WRITE(iutarl,'(20x,f6.0,i3,20(1x,a4))')  &  
           1.0 ,n0,(vchar0(n),n=1,n0)

! upper level information
  DO k = kx, (kx-ka+1), -1
     WRITE(iutarl,'(20x,f6.4,i3,20(1x,a4))')  &
           arlsig(k),n1,(vchar1(n),n=1,n1)
  END DO

  CLOSE (iutarl)
  RETURN

END SUBROUTINE setgrid  
