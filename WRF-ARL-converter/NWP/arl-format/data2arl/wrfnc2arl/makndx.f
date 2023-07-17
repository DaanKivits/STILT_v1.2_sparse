! Use emacs f90-mode:   -*- mode:f90  -*-
module makendx_module

! $Id: makndx.f,v 1.3 2013/07/31 13:58:24 trn Exp $

contains    

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAKNDX           CREATE THE INDEX RECORD CONFIGURATION FILE
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
!     UNIT iucfg - one of CFG_(MASS|STAX|STAY).WRF, defines the output file grid and structure
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

!USAGE: 
      SUBROUTINE MAKNDX(handle, igtyp, nlevtyp, p_sfc, p_half, p_full, &    !
              & ngridtyp, p_mass, p_stax, p_stay, &                       !
              & vnames, nvars, nsigl, ndim1l, ndim1v, &                   !
              & sigl, iprint, iucfg, ig, nxp, nyp, &                      !
              & nvar, lvl, nlvl, combine_grids, combine_levs, mdlid, geoloc_type)             !

      use setmap_module

      implicit none

      include 'netcdf.inc'

!   INPUT ARGUMENT LIST:

      INTEGER, intent(in) :: handle, igtyp, nlevtyp, p_sfc, p_half, p_full, &
           & ngridtyp, p_mass, p_stax, p_stay, ndim1l, ndim1v
!     handle - netcdf file unit number
!     igtyp - grid type (one of the ngridtyp values: p_mass, p_stax, p_stay)
!     nlevtyp - number of level types: p_sfc, p_half, p_full
!     ndim1l, ndim1v - leading dimensions (max values) for levels and variables

      character (len=4), intent(in) :: vnames(ndim1v,nlevtyp,ngridtyp), mdlid, geoloc_type
!     Vnames - character array of variable names
!     mdlid - model ID for output in cfg file
!     geoloc_type - which mapping routines to use

      integer, intent(in) :: nvars(nlevtyp,ngridtyp), nsigl(nlevtyp), &
           & sigl(ndim1l,nlevtyp), iprint, iucfg, ig(ngridtyp)
!     sigl - level values (including sfc), by level type
!     iprint - printout control
!     iucfg - unit number to use for configuration file (reusable)
!     IG - grid identification number

      logical, intent(in) :: combine_grids, combine_levs
!     flags controlling whether to combine grids, levels, and whether to create index

!   OUTPUT ARGUMENT LIST (output for current igtyp, unchanged for others):
      integer, intent(inout) :: nxp(ngridtyp), nyp(ngridtyp)
!     NXP,NYP - output grid dimensions

      integer, intent(inout) :: nvar(2*ndim1l,ngridtyp), lvl(2*ndim1l,ngridtyp), nlvl(ngridtyp)
!     NLVL - number of data levels in output file, by grid type
!     nvar - number of variables in output file, by level and grid type
!     lvl - height of each output level, by grid type
!           (merges sfc, full, and half levels as needed)

!   Local data:
      LOGICAL :: FTEST
      REAL   :: grids(12), sigma
      integer :: i, ii, ii1, ii2, iladd, ilev, kl, kladd, lev2, nl, nv
      integer :: nxp_out, nyp_out, igtyp_out, ierr

!     arrays to hold variable selection information

      CHARACTER (len=4) :: VCHAR0(2*ndim1v), VCHAR1(2*ndim1v,2*ndim1l), fname_grid
!     VCHAR0 - character array of surface field identification
!     VCHAR1 - character array of upper level field identifications, by level

      character (len=80) :: fname

      call setmap(iprint .ge. 8,handle,nxp(igtyp),nyp(igtyp),grids,ierr,geoloc_type)
      if (ierr .ne. 0) stop 'Error in setmap'
!     constrain grid point 1,1 longitude to [-180, +180] (needed ??)      
!       if (grids(11) .gt.  180.0) grids(11) = grids(11) - 360.0
!       if (grids(11) .lt. -180.0) grids(11) = grids(11) + 360.0

!     write the packer configuration file
      
      igtyp_out=igtyp
      if (combine_grids) igtyp_out = p_mass
      if (igtyp_out .eq. p_mass) fname_grid='MASS'
      if (igtyp_out .eq. p_stax) fname_grid='STAX'
      if (igtyp_out .eq. p_stay) fname_grid='STAY'
      fname = 'CFG_' // fname_grid // '.WRF'

      if (iprint .ge. 2) write (*,'(a,i10,1x,a)') &
           & 'Writing cfg info to unit, file=',iucfg,trim(fname)
      OPEN(iucfg,file=fname)
      WRITE(iucfg,'(20X,A4)') mdlid

!     grid number 99 and 1 for sigma coordinate system
      WRITE(iucfg,'(20X,I4)') IG(igtyp_out), 1

      WRITE(iucfg,'(20X,F10.2)')(GRIDS(I),I=1,12)

! level and variable information
      if (combine_grids) then
         ii1=1
         ii2=ngridtyp
      else
         ii1=igtyp
         ii2=igtyp
      endif
      nvar(:,igtyp_out) = 0
      nlvl = 0                  !Added by DVM and DWen
! sfc:
! sfc variables, if any:
      do ii=ii1,ii2
         do i=1,nvars(p_sfc,ii)
            vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_sfc,ii)
         end do
         nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_sfc,ii)
! add lowest full level variables, if any
         if (nvars(p_full,igtyp) .gt. 0) then
            do i=1,nvars(p_full,ii)
               vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_full,ii)
            end do
            nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_full,ii)
         end if
      end do
! if neither sfc nor full variables there, use lowest half-level
! (will be missing-filled)
      if (nvar(1,igtyp_out) .eq. 0) then
         do ii=ii1,ii2
            do i=1,nvars(p_half,ii)
               vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_half,ii)
            end do
            nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_half,ii)
         enddo
      endif
      if (iprint .ge. 4) write (*,'(a,(10(a4,1x)))') 'vchar0=',(vchar0(i),i=1,nvar(1,igtyp_out))

! continue here

! First level is at surface:
      kl = 1
      lvl(kl,igtyp_out) = 10000
! Upper levels:
      lev2 = max(nsigl(p_half),nsigl(p_full))
      do ilev=1,lev2-1
         kladd=1
! For half-levels, loop over all grid types that are to be combined:
         do ii=ii1,ii2
            if (nvars(p_half,ii) .gt. 0) then
! Increment output level counter
               kl = kl + kladd
! Half-level variables:
               do i=1,nvars(p_half,ii)
                  vchar1(nvar(kl,igtyp_out)+i,kl) = vnames(i,p_half,ii)
               end do
               nvar(kl,igtyp_out) = nvar(kl,igtyp_out) + nvars(p_half,ii)
! Disable level incrementing for subsequent grid and level types
               kladd = 0
               if (.not. combine_levs) then
                  !           Keep half-levels as separate levels:
                  iladd=1
                  lvl(kl,igtyp_out) = sigl(ilev+iladd,p_half)
                  if (iprint .ge. 4) write (*,'(a,4i10)') 'half ilev,kl,lvl,nvar=', &
                       & ilev,kl,lvl(kl,igtyp_out),nvar(kl,igtyp_out)
                  if (lvl(kl,igtyp_out) .ge. lvl(kl-1,igtyp_out)) stop 'levels out of order'
               end if
            end if
         end do

! Enable level incrementing for full levels if not combining levels
         if (.not. combine_levs) kladd = 1
! For full-levels, loop over all grid types that are to be combined:
         do ii=ii1,ii2
            if (nvars(p_full,ii) .gt. 0) then
! Increment output level counter
               kl = kl + kladd
! Half-level variables:
               do i=1,nvars(p_full,ii)
                  vchar1(nvar(kl,igtyp_out)+i,kl) = vnames(i,p_full,ii)
               end do
               nvar(kl,igtyp_out) = nvar(kl,igtyp_out) + nvars(p_full,ii)
! Disable level incrementing for subsequent grid and level types
               kladd = 0
               iladd=1
               lvl(kl,igtyp_out) = sigl(ilev+iladd,p_full)
               if (iprint .ge. 4) &
                    & write (*,'(a,4i10)') 'full ilev,kl,lvl,nvar=', &
                    & ilev,kl,lvl(kl,igtyp_out),nvar(kl,igtyp_out)
               if (lvl(kl,igtyp_out) .ge. lvl(kl-1,igtyp_out)) stop 'levels out of order'
            end if
         end do

      end do
       
      nlvl(igtyp_out) = kl

      if (combine_grids) then
         if (igtyp .ne. p_mass) stop 'makndx: igtyp internal logic error'
         nxp_out=nxp(igtyp)+1
         nyp_out=nyp(igtyp)+1
      else
         nxp_out=nxp(igtyp)
         nyp_out=nyp(igtyp)
      endif
      WRITE(iucfg,'(20X,I4)')nxp_out,nyp_out,nlvl(igtyp_out)

      DO NL=1,NLVL(igtyp_out)
         sigma=lvl(nl,igtyp_out)/10000.0
         if(nl.eq.1)then
!     sfc level information
            WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
            sigma,nvar(nl,igtyp_out),(vchar0(nv),nv=1,nvar(nl,igtyp_out))
         else
!     upper level information
            WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
            sigma,nvar(nl,igtyp_out),(vchar1(nv,nl),nv=1,nvar(nl,igtyp_out))
         end if
      END DO
      CLOSE (iucfg)
      RETURN
      end SUBROUTINE MAKNDX
    end module makendx_module
    
