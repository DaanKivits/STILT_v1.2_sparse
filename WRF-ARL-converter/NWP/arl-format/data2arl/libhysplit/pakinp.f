!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKINP           PAcK INPut converts char*1 to real
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PACK INPUT DOES THE CONVERSION OF CHAR*1 PACKED ARRAY TO
!   A REAL*4 DATA ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISION: 14 Feb 1997 (RRD)
!                  04 Feb 2000 (RRD) - Sun F90 compatibility option
!                  29 Sep 2000 (RRD) - fortran90 upgrade
!                  02 Apr 2001 (RRD) - added function definition
!                  09 Sep 2002 (RRD) - fortran coding standards
!                  25 Feb 2003 (RRD) - ichar generalization
!
! USAGE:  CALL PAKINP(RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
!
!   INPUT ARGUMENT LIST:      see below
!   OUTPUT ARGUMENT LIST:     see below
!   INPUT FILES:              none
!   OUTPUT FILES:             none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PAKINP(RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  REAL,          INTENT(OUT)   :: rvar (:,:)     ! real data unpacked
  CHARACTER(1),  INTENT(IN)    :: cvar (:)       ! packed input of NX*NY
  INTEGER,       INTENT(IN)    :: nx,ny          ! size of input array  
  INTEGER,       INTENT(IN)    :: nx1,ny1        ! optional sub-grid left edge 
  INTEGER,       INTENT(IN)    :: lx,ly          ! length of sub-grid
  REAL,          INTENT(IN)    :: prec           ! precision of packed data 
  INTEGER,       INTENT(IN)    :: nexp           ! packing scaling exponent
  REAL,          INTENT(IN)    :: var1           ! value of array at (1,1)
  INTEGER,       INTENT(INOUT) :: ksum           ! rotating checksum 

!-------------------------------------------------------------------------------
! internal variable definitions
!-------------------------------------------------------------------------------

  REAL                         :: scexp,rold,rnew
  INTEGER                      :: i,ii,j,jj,k

!-------------------------------------------------------------------------------
! Some Fortran compilers the ICHAR function is undefined for values greater
! Generic Fortran and C code functions for conversions can be used (ncep mova2i)
! In Fortran:  JCHAR=ISHFT(TRANSFER(MYCHR,JCHAR),8-BIT_SIZE(JCHAR))
! In C:        INT MOVA2I(UNSIGNED CHAR *A) { RETURN (INT) (*A); }
! If ICHAR always returns a signed integer, then ICHAR can be replaced with:
!-------------------------------------------------------------------------------
  CHARACTER(1)  :: mychr    
  INTEGER       :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
!-------------------------------------------------------------------------------

  SCEXP= 1.0 / 2.0**(7-NEXP)  ! scaling exponent

!-------------------------------------------------------------------------------
! unpack initial value for each row of column 1
!-------------------------------------------------------------------------------

  ROLD=VAR1
  DO J=1,NY
     K=(J-1)*NX+1                                 ! position at column 1
     RNEW=(FLOAT(JCHAR(CVAR(K))-127)*SCEXP)+ROLD  ! value from prevous row
     ROLD=RNEW
     JJ=J-NY1+1                                   ! index in output sub-grid
     IF(JJ.GE.1.AND.JJ.LE.LY)RVAR(1,JJ)=RNEW      ! sub-grid at left edge
  END DO

!-------------------------------------------------------------------------------
! only unpack within J-subgrid
!-------------------------------------------------------------------------------

  DO J=NY1,(NY1+LY-1)
     JJ=J-NY1+1                             ! sub-grid array (1 to LY)

     ROLD=RVAR(1,JJ)
     DO I=2,(NX1+LX-1)                      ! unpack I from 1 to I sub-grid max
        K=(J-1)*NX+I
        RNEW=(FLOAT(JCHAR(CVAR(K))-127)*SCEXP)+ROLD
        ROLD=RNEW
        II=I-NX1+1                          ! sub-grid array element (1 to LX)
        IF(ABS(RNEW).LT.PREC)RNEW=0.0       ! check precision for true zero
        IF(II.GE.1.AND.II.LE.LX)RVAR(II,JJ)=RNEW
     END DO
  END DO

!-------------------------------------------------------------------------------
! only do full-grid checksum when KSUM=0
!-------------------------------------------------------------------------------

  IF(KSUM.NE.0)RETURN

  DO K=1,(NX*NY)
     KSUM=KSUM+JCHAR(CVAR(K))
     IF(KSUM.GE.256)KSUM=KSUM-255    ! sum carries over the eighth bit add one
  END DO

END SUBROUTINE pakinp
