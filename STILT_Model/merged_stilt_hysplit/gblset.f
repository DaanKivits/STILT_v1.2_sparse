!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  GBLSET           GloBaL grid SETup   
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:01-09-03
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   SETUPS THE PARAMETERS FOR A GLOBAL LATITUDE LONGITUDE GRID
!   THAT MAY BE USED FOR COMPUTATIONS INSTEAD OF THE CONFORMAL   
!   PROJECTION SYSTEM.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Mar 2001 (RRD) - initial version
!                 29 Aug 2001 (RRD) - simultaneous multiple meteo
!                 22 Feb 2002 (RRD) - added maxgrid initialization
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 10 Aug 2006 (RRD) - prime meridian test
!
! USAGE:  CALL GBLSET(KG,KT)  
! 
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none
!   OUTPUT FILES:            none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE GBLSET(KG,KT)

  use module_defgrid ! meteorology grid and file
  IMPLICIT NONE
 

  INTEGER, INTENT(IN) :: kg    ! active grid number    
  INTEGER, INTENT(IN) :: kt    ! active time number    

  REAL    :: DSX,DSY,CLAT1,CLON1,CLAT2,CLON2,DLON 


  GRID(KG,KT)%LATLON=.FALSE.
 ! GRID(KG,KT)%GLOBAL=.FALSE.
  GRID(KG,KT)%GLOBAL='no' ! default
!  GRID(KG,KT)%GBLDAT=.FALSE.
   GRID(KG,KT)%GBLDAT="no"
  GRID(KG,KT)%PRIME =.FALSE.

! determine if this is a lat-lon grid
  IF(GRID(KG,KT)%SIZE.EQ.0.0)THEN
     GRID(KG,KT)%LATLON=.TRUE.
     CALL GBLDXY(KG,KT,1.0,1.0,DSX,DSY)
     GRID(KG,KT)%SIZE=DSY 
  ELSE
     RETURN
  END IF

! find the corner points
  CLAT1=GRID(KG,KT)%SYNC_LAT
  CLON1=GRID(KG,KT)%SYNC_LON
  CLAT2=GRID(KG,KT)%POLE_LAT
  CLON2=GRID(KG,KT)%POLE_LON

! grid spacing
  DLON=GRID(KG,KT)%REF_LON
  write(*,*)CLAT1,CLAT2,CLON1,CLON2,DLON

! determine if the grid is global
  IF( ( (CLON2+DLON-CLON1.EQ.360.0).OR.   &
     (CLON2+DLON-CLON1.EQ.  0.0) ) .AND.  &
      CLAT2-CLAT1.EQ.180.0) THEN
   !   GRID(KG,KT)%GBLDAT=.TRUE.  
       GRID(KG,KT)%GBLDAT="gl" 

 !tk(20160317)
  ELSE IF( ( (CLON2+DLON-CLON1.EQ.360.0).OR.   &
     (CLON2+DLON-CLON1.EQ.  0.0) ) .AND.  &
      CLAT2+DLON.EQ.90.5) THEN
   !   GRID(KG,KT)%GBLDAT=.TRUE.  
       GRID(KG,KT)%GBLDAT="nh"
 !tk(20160317)
  ELSE IF( ( (CLON2+DLON-CLON1.EQ.360.0).OR.   &
     (CLON2+DLON-CLON1.EQ.  0.0) ) .AND.  &
      CLAT1.EQ. -90.0) THEN
   !   GRID(KG,KT)%GBLDAT=.TRUE.  
       GRID(KG,KT)%GBLDAT="sh"
  END IF

! determine if a non-global grid is about the prime 
  IF(CLON2.GE.0.0.AND.CLON1.LT.0.0) GRID(KG,KT)%PRIME=.TRUE.  

END SUBROUTINE gblset

