! Use emacs f90-mode:   -*- mode:f90  -*-
module datcnv_module

! $Id: datcnv.f,v 1.2 2007/11/20 22:37:49 trn Exp $

contains    

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DATCNV           CONVERT UNITS FOR ALL ELEMENTS OF THE DATA
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            IS USED TO CONVERT EACH ELEMENT OF THE OUTPUT DATA GRID
!            TO CONFORM WITH UNIT CONVENTIONS USED IN OTHER APPLICATIONS
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 - RRD
!
! USAGE:  CALL DATCNV(RVAR,NXP,NYP,CNVRT)
!   INPUT ARGUMENT LIST:
!     RVAR - real data array
!     NXP,NYP - dimensions of the array
!     CNVRT - conversion factor
!   OUTPUT ARGUMENT LIST:
!     RVAR - real data array after conversion
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE datcnv(rvar,nxp,nyp,nlev,ntimes,itime,cnvrt,add)

      implicit none

      integer, intent(in) :: nxp, nyp, nlev, ntimes, itime
      real, intent(in) :: cnvrt, add
      REAL, intent(inout) ::   RVAR(NXP,NYP,nlev,ntimes)

      integer :: i,j,k

      do k=1,nlev
      DO J=1,NYP
      DO I=1,NXP

         rvar(i,j,k,itime)=(rvar(i,j,k,itime) + add)*cnvrt

      END DO
      END DO
      END DO
      RETURN
      end SUBROUTINE datcnv
    end module datcnv_module
    
