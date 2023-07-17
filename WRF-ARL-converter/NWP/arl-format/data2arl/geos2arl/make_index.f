! Use emacs f90-mode:   -*- mode:f90  -*-
module make_index_module

  implicit none
  private
  public :: make_index

  ! $Id: make_index.f,v 1.2 2013/04/01 16:45:26 trn Exp $

contains    

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  MAKE_INDEX           CREATE THE INDEX RECORD CONFIGURATION FILE
  !
  ! Written by Thomas Nehrkorn, AER, Inc.
  ! Based on:
  !   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
  !
  ! ABSTRACT:  CREATES THE CONFIGURATION FILE FOR THE OUTPUT DATA WHICH
  !            DEFINES THE GRID SYSTEM AS WELL AS ALL VARIABLES THAT ARE
  !            TO BE WRITTEN TO THE OUTPUT FILE.
  !
  ! USAGE:  see below
  !   INPUT FILES:
  !     NONE
  !   OUTPUT FILES:
  !     UNIT iucfg - one of CFG_timestamp.WRF, defines the output file grid and structure
  !
  ! ATTRIBUTES:
  !   LANGUAGE: FORTRAN 90
  !   MACHINE:  portable
  !
  !$$$

  !USAGE: 
  SUBROUTINE MAKE_INDEX(nlevtyp, p_sfc, p_lay, &
       vnames, nvars, nsigl, ndim1l, ndim1v, &
       sigl, iprint, iucfg, nx, ny, mdlid, grids, &
       nvar, lvl, nlvl)

    implicit none

    !   INPUT ARGUMENT LIST:

    INTEGER, intent(in) :: nlevtyp, p_sfc, p_lay, ndim1l, ndim1v
    !     nlevtyp - number of level types: p_sfc, p_lay
    !     ndim1l, ndim1v - leading dimensions (max values) for levels and variables

    character (len=4), intent(in) :: vnames(ndim1v,nlevtyp), mdlid
    !     Vnames - character array of variable names
    !     mdlid - model ID for output in cfg file

    integer, intent(in) :: nvars(nlevtyp), nsigl(nlevtyp), sigl(ndim1l,nlevtyp), iprint, iucfg
    !     sigl - level values (including sfc), by level type
    !     iprint - printout control
    !     iucfg - unit number to use for configuration file (reusable)

    integer, intent(in) :: nx, ny
    !     NX,NY - output grid dimensions

    REAL, intent(in)   :: grids(12)

    !   OUTPUT ARGUMENT LIST:
    integer, intent(inout) :: nvar(2*ndim1l), lvl(2*ndim1l), nlvl
    !     NLVL - number of data levels in output file
    !     nvar - number of variables in output file, by level
    !     lvl - height of each output level
    !           (merges sfc, and lay levels as needed)

    !   Local data:
    integer :: i, ilev, kl, kladd, lev2, nl, nv
    real :: sigma

    !     arrays to hold variable selection information

    CHARACTER (len=4) :: VCHAR0(2*ndim1v), VCHAR1(2*ndim1v,2*ndim1l)
    !     VCHAR0 - character array of surface field identification
    !     VCHAR1 - character array of upper level field identifications, by level

    !     write the packer configuration file

    WRITE(iucfg,'(20X,A4)') mdlid

    !     grid number 99 and 1 for sigma coordinate system
    WRITE(iucfg,'(20X,I4)') 99, 1

    !     mapping info
    WRITE(iucfg,'(20X,F10.2)') (GRIDS(I),I=1,12)

    ! level and variable information
    nvar(:) = 0
    ! sfc:
    ! sfc variables, if any:
    do i=1,nvars(p_sfc)
       vchar0(nvar(1)+i) = vnames(i,p_sfc)
    end do
    nvar(1) = nvar(1) + nvars(p_sfc)
    ! if no sfc variables there, use lowest lay-level
    ! (will be missing-filled)
    if (nvar(1) .eq. 0) then
       do i=1,nvars(p_lay)
          vchar0(nvar(1)+i) = vnames(i,p_lay)
       end do
       nvar(1) = nvar(1) + nvars(p_lay)
    endif
    if (iprint .ge. 4) write (*,'(a,(10(a4,1x)))') 'vchar0=',(vchar0(i),i=1,nvar(1))

    ! First level is at surface:
    kl = 1
    lvl(kl) = 10000
    ! Upper levels:
    lev2 = nsigl(p_lay)
    do ilev=2,lev2
       kladd=1
       ! For layers:
       if (nvars(p_lay) .gt. 0) then
          ! Increment output level counter
          kl = kl + kladd
          ! Lay-level variables:
          do i=1,nvars(p_lay)
             vchar1(nvar(kl)+i,kl) = vnames(i,p_lay)
          end do
          nvar(kl) = nvar(kl) + nvars(p_lay)
          lvl(kl) = sigl(ilev,p_lay)
          if (iprint .ge. 4) write (*,'(a,4i10)') 'lay ilev,kl,lvl,nvar=', &
               & ilev,kl,lvl(kl),nvar(kl)
          if (lvl(kl) .ge. lvl(kl-1)) write (*,*) 'make_index warning: levels out of order ?'
       end if
    end do

    nlvl = kl

    WRITE(iucfg,'(20X,I4)')nx,ny,nlvl

    DO NL=1,NLVL
       sigma=lvl(nl)/10000.0
       if(nl.eq.1)then
          !     sfc level information
          WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
               sigma,nvar(nl),(vchar0(nv),nv=1,nvar(nl))
       else
          !     upper level information
          WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
               sigma,nvar(nl),(vchar1(nv,nl),nv=1,nvar(nl))
       end if
    END DO
    RETURN
  end SUBROUTINE MAKE_INDEX
end module make_index_module
