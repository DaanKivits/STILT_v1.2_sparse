! Use emacs f90-mode:   -*- mode:f90  -*-
module xtract_module

! $Id: xtract.f,v 1.9 2016/04/18 14:48:55 trn Exp $

contains
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EXTRACT          MAIN ROUTINE FOR DECODING WRF GRIB DATA
! Written by Thomas Nehrkorn, AER, Inc.
!   Based on:
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  EXTRACTS all RECORDs FROM INPUT netcdf FILE that correspond to
!            the specified variables, levels, and writes them
!            TO THE OUTPUT FILE
!
! USAGE:  see below
!   INPUT FILES:
!     HANDLE defines input data file
!   OUTPUT FILES:
!     unit 20(iumass) DATA_MASS.WRF - ARL packed data output file for mass points
!     unit 21(iustax) DATA_STAX.WRF - ARL packed data output file for x-staggered points
!     unit 22(iustay) DATA_STAY.WRF - ARL packed data output file for y-staggered points
!     unit 30(iucfg) - used for CFG_MASS.WRF, CFG_STAX.WRF, CFG_STAY.WRF (temp file area)
!     unit 40(iutime) WRFTIME - time indicator file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

      SUBROUTINE xtract(handle, iumass, iustax, iustay, iucfg, iutime, iprint, &
                 & nlevtyp, p_sfc, p_half, p_full, ngridtyp, p_mass, p_stax, p_stay, &
                 & nvars, ndim1v, &
                 & vlabels, lchar, vnames, vcnvrt, vadd,&
                 & combine_grids, combine_levs, mdlid, use_16, itime, geoloc_type )

      use makendx_module
      use getdim_module
      use gettime_module
      use setvar_module
      use datini_module
      use datcnv_module
      use copy_arr_module
      use hybrid_vertical_coords_module

      implicit none 

      include 'netcdf.inc'
      include 'pakset_16.INC'

!   INPUT ARGUMENT LIST:
      INTEGER, intent(in) :: handle, lchar
!     HANDLE - defines input file to directIO routines
      integer, intent(in) :: iumass, iustax, iustay, iucfg, iutime, iprint !units, printout ctl
      integer, intent(in) :: nlevtyp, p_sfc, p_half, p_full !level type dimension and values
      integer, intent(in) :: ngridtyp, p_mass, p_stax, p_stay !grid type dimension and values
      integer, intent(in) :: nvars(nlevtyp,ngridtyp)
!     nvars - active number of array elements of vnames, by level and grid type
      integer, intent(in) :: ndim1v !leading dimensions (max values) for variables
      character (len=lchar), intent(in) :: vlabels(ndim1v,nlevtyp,ngridtyp) !netcdf variable names
      real, intent(in) :: vcnvrt(ndim1v,nlevtyp,ngridtyp) !conversion factor
      real, intent(in) :: vadd(ndim1v,nlevtyp,ngridtyp) !add to input before conversion
      character (len=4), intent(in) :: vnames(ndim1v,nlevtyp,ngridtyp) !variable names
      logical, intent(in) :: combine_grids, combine_levs, use_16
      integer, intent(in) :: itime !time period numb er to extract
      integer :: nxy_factor, nxy_arg, itime_local
      character (len=4), intent(inout) :: mdlid
      character (len=4), intent(in) :: geoloc_type

!     number of levels, variables, and array limits
      integer :: ndim1l !leading dimension (max values) for levels 
      integer, allocatable :: nsigl(:) ! active number of levels, by level type
      integer, allocatable :: sigl(:,:) !level values (including sfc), by level type
      character (len=lchar) :: vlabel, label !netcdf variable names
      real,allocatable :: xsigl(:,:),soildepths(:,:)
      real ::  sig_factor, min_depth
      integer :: iladd_save, iladd_soil, nsoil
      integer :: ndim,dimlen(nf_max_var_dims),n3d,varid
      integer :: MAXX, MAXY, MAXB
      integer :: MAXC

!     time decoding
      character (len=14) :: init_dattim,current_dattim
      integer :: iyyyy, lead_time, lead_units, iss, ierr_tim, any_ierr_tim

!     remap array from one to two dimensions
      integer :: nxp(ngridtyp), nyp(ngridtyp), nxy(ngridtyp)
      integer :: nxp_pakset, nyp_pakset, nzp_pakset
      REAL, allocatable :: SVAR(:,:)
!     unpacked output array
      REAL, allocatable :: rvar_mass(:,:,:,:),rvar_stax(:,:,:,:),rvar_stay(:,:,:,:)
      REAL, allocatable :: rvar_out(:,:)
      INTEGER :: shape(2,ngridtyp),nlvl(ngridtyp)
      INTEGER , allocatable :: nvar(:,:), lvl(:,:)

!     packed output array
      CHARACTER (len=1), allocatable ::  CVAR(:)

!     default information (grid id, output records)
      character (len=80) :: fname
      character (len=4) :: fname_grid
      integer :: iudata
      integer :: ig(ngridtyp)
      integer :: krec(ngridtyp), krec_out(ngridtyp)
      integer, save :: icall=0
      integer :: status, dimid

      logical :: have_index(ngridtyp), have_init(ngridtyp)
      integer :: i
      integer :: ierr, ntimes
      integer :: kbyte
      integer :: klen
      integer :: koff
      integer :: kvarb
      integer :: ktype
      integer :: level
      logical :: lmatch
      integer :: igtyp, igtyp_out
      integer :: iltyp
      integer :: ivar
      integer :: ilev
      integer :: iladd
      integer :: kret
      integer :: iyr
      integer :: imo
      integer :: ida
      integer :: ihr
      integer :: imn
      integer :: ifh
      integer :: kl
      integer :: idatlin, ndatlin, nlev_in, nlev_var


!-------------------------------------------------------------------------------
! only required when dealing with some F90 compiler
! replace ICHAR below with internally defined JCHAR function
  CHARACTER(1)                 :: mychr    
  INTEGER                      :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
!-------------------------------------------------------------------

      if (icall .eq. 0) then
         ig(:)=99
         krec(:)=0
         krec_out(:)=0
         have_index(:) = .FALSE.
         have_init(:) = .FALSE.
      else
         stop 'multiple netcdf files not supported, call it one at a time'
      endif
      icall = icall + 1

!#    generate valid date/time:
      itime_local=itime
      call gettime(iprint .ge. 8,handle,itime_local,iyyyy,imo,ida,ihr,imn,iss,iyr,ntimes,ierr)
      if (ierr .ne. 0) stop 'Error getting time'
      if (iprint .ge. 2) write (*,'(a,i4,a,i4)') &
           & 'Using time period',itime,' out of a total possible',ntimes
      if (iprint .ge. 2) write (*,'(a,i4,5(a1,i2.2))') &
           & 'Decoded time: ',iyyyy,'-',imo,'-',ida,'_',ihr,':',imn,':',iss
      ifh = 0

! Get level information
      allocate(nsigl(nlevtyp))
      nsigl(:) = 0
      nsigl(p_sfc) = 1
      ndim1l=0
      sig_factor=-1.
      do i=1,2
         if (i .eq. 1) then
            iltyp=p_full
            iladd=0
            label = 'bottom_top_stag'
            vlabel = 'ZNW'
         else
            iltyp=p_half
            iladd=1
            label = 'bottom_top'
            vlabel = 'ZNU'
         endif
         call getdim(iprint .ge. 8,handle,label,ndatlin,ierr)
         if (ierr .ne. 0) stop 'Error getting dim for sigma/eta'
         if (ndim1l .eq. 0) then
            ndim1l=ndatlin
            if (iprint .ge. 2) write (*,*) 'Allocating sigl with ndatlin= ',ndatlin
            allocate(xsigl(ndatlin,ntimes), &
                 & sigl(ndim1l,nlevtyp), nvar(2*ndim1l,ngridtyp), lvl(2*ndim1l,ngridtyp), stat=ierr)
            if (ierr .ne. 0) stop 'Allocate sigl...'
            sigl(1,p_sfc)=10000
            sigl(1,p_half)=10000
         endif
         nsigl(iltyp) = ndatlin+iladd
         call setvar(iprint .ge. 8,handle,vlabel,varid,n3d,ndim,dimlen,ierr)
         if (ierr .ne. 0) stop 'setvar error for sigma/eta values'
         if (ndim .ne. 2 .or. dimlen(1) .ne. ndatlin) stop 'inconsistent eta dimensions'
         status = nf_get_var_real(handle,varid,xsigl)
         if (status .ne. nf_noerr) stop 'error getting sigma/eta values'
         if (sig_factor .le. 0) then
            if (xsigl(1,1) .lt. 100) then
               sig_factor=10000.
            else
               sig_factor=1.
            endif
         endif
         do idatlin=1,ndatlin
            sigl(idatlin+iladd,iltyp)=nint(sig_factor*xsigl(idatlin,1))
         end do
      end do
      
! check that full/half levels are provided if needed, and that they are consistent
      if (nsigl(p_half) .ne. nsigl(p_full)) &
           & stop 'inconsistent half and full levels specified'

      if (iprint .ge. 3) then
         if (nsigl(p_half) .gt. 1) write (*,'(a/,(2i10))') 'Specified half-levels:', &
              & (i,sigl(i,p_half),i=1,nsigl(p_half))
         if (nsigl(p_full) .gt. 1) write (*,'(a/,(2i10))') 'Specified full-levels:', &
              & (i,sigl(i,p_full),i=1,nsigl(p_full))
         do igtyp=1,ngridtyp
            do iltyp=1,nlevtyp
               write (*,'(a,3i10)') 'Grid type, level type, nvars=', &
                    & igtyp,iltyp,nvars(iltyp,igtyp)
               if (nvars(iltyp,igtyp) .gt. 0) write (*,'(a4,1x,a4,1x,a40,2a20)') &
                    & 'ivar','name',' netcdf label ',' conversion factor ',' add before conv '
               do i=1,nvars(iltyp,igtyp)
                  write (*,'(i4,1x,a4,1x,a40,2g20.5)') i,vnames(i,iltyp,igtyp), &
                       & trim(vlabels(i,iltyp,igtyp)),vcnvrt(i,iltyp,igtyp),vadd(i,iltyp,igtyp)
               end do
            end do
         end do
      end if

! processing for hybrid vertical coordinate
      call hybrid_vertical_coords(iprint, handle, nsigl(p_full), ntimes, iucfg, mdlid, ierr)
      if (ierr .ne. 0) stop 'Error in hybrid_vertical_coords'
      
! get array dimensions
      igtyp = p_mass
      label = 'west_east'
      call getdim(iprint .ge. 8,handle,label,nxp(igtyp),ierr)
      if (ierr .ne. 0) stop 'Error getting west_east dim'
      nxp(p_stay) = nxp(p_mass)
      label = 'south_north'
      call getdim(iprint .ge. 8,handle,label,nyp(igtyp),ierr)
      if (ierr .ne. 0) stop 'Error getting south_north dim'
      nyp(p_stax) = nyp(p_mass)
      
      igtyp = p_stax
      label = 'west_east_stag'
      call getdim(iprint .ge. 8,handle,label,nxp(igtyp),ierr)
      if (ierr .ne. 0) stop 'Error getting west_east_stag dim'

      igtyp = p_stay
      label = 'south_north_stag'
      call getdim(iprint .ge. 8,handle,label,nyp(igtyp),ierr)
      if (ierr .ne. 0) stop 'Error getting south_north_stag dim'

      label = 'soil_layers_stag'
      call getdim(iprint .ge. 8,handle,label,nsoil,ierr)
      if (ierr .ne. 0) stop 'Error getting soil_layers_stag'

      vlabel = 'ZS'
      call setvar(iprint .ge. 8,handle,vlabel,varid,n3d,ndim,dimlen,ierr)
      if (ierr .ne. 0) then
         write (*,*) 'WARNING: setvar error for soil depths: ZS, not checking soil dims'
      else
         if (ndim .ne. 2 .or. dimlen(1) .ne. nsoil .or. dimlen(ndim) .ne. ntimes) then
            write (*,*) 'Unexpected dimensions for soil layers: ndim= ', &
                 & ndim,' dimlen =',dimlen(1:ndim)
            stop 'Unexpected dimensions for soil layers'
         end if
      end if
      allocate(soildepths(nsoil,ntimes), stat=ierr)
      if (ierr .ne. 0) stop 'Error allocating soildepths'
      status = nf_get_var_real(handle,varid,soildepths)
      if (status .ne. nf_noerr) then
         write (*,*) 'WARNING: Error getting soildepths (ZS) from netcdf file.'
         write (*,*) 'Assuming they are increasing in absolute value, make sure this is correct'
         do i=1,nsoil
            soildepths(i,:) = i
         end do
      endif
      iladd_soil = 0
      min_depth = 1e10
      do i=1,nsoil
         if (abs(soildepths(i,1)) .lt. min_depth) then
            iladd_soil = i-1
            min_depth = abs(soildepths(i,1))
         end if
      end do
      
      if (iprint .ge. 4) write (*,*) 'nxp= ',nxp,' nyp= ',nyp,' nsoil=',nsoil,' iladd_soil=',iladd_soil
      maxx = nxp(p_stax)
      maxy = nyp(p_stay)

!     Initial array sizes
      nxy_factor = 1
      if (use_16) nxy_factor = 2
      MAXC=nxy_factor * maxx * maxy

      deallocate (svar, cvar, stat=ierr)
      allocate (cvar(maxc), stat=ierr)
      if (ierr .ne. 0) stop 'Initial allocate buff...'

! loop over igtyp, iltyp, ivar, ilev
!     check if variable present in selection table
      IGTYP_LOOP: do igtyp=1,ngridtyp
         igtyp_out=igtyp
         if (combine_grids) then
            igtyp_out = p_mass
         else
            write (*,*) 'combine_grids=.FALSE. is not supported'
            stop 'combine_grids=.FALSE. is not supported'
!       TBD would require changes in makndx (map info for staggered grids), and calls to makndx
!           for each igtyp
         end if
         
         if (igtyp_out .eq. p_mass) then
            iudata=iumass
         elseif (igtyp_out .eq. p_stax) then
            iudata=iustax
         elseif (igtyp_out .eq. p_stay) then
            iudata=iustay
         endif
!     create an index record for pakrec and perform initializations: only supporting igtyp_out=pmass at the moment
         if (igtyp .eq. 1) &
              & call makndx(handle, igtyp, nlevtyp, p_sfc, p_half, p_full, &
              & ngridtyp, p_mass, p_stax, p_stay, &
              & vnames, nvars, nsigl, ndim1l, ndim1v, &
              & sigl, iprint, iucfg, ig, nxp, nyp, &
              & nvar, lvl, nlvl, combine_grids, combine_levs, mdlid, geoloc_type)

!     perform initializations for this grid typ
         if (.not. have_init(igtyp)) then
            shape(1,igtyp)=nxp(igtyp)
            shape(2,igtyp)=nyp(igtyp)
            nxy(igtyp) = nxp(igtyp)*nyp(igtyp)
            if (combine_grids .and. .not. have_init(igtyp_out)) then
               if (igtyp .ne. igtyp_out) stop 'xtract: igtyp internal logic error'
               nxy(igtyp) = (nxp(igtyp)+1)*(nyp(igtyp)+1)
               allocate(rvar_out(nxp(igtyp)+1,nyp(igtyp)+1),stat=ierr)
               if (ierr .ne. 0) stop 'Allocate rvar_out...'
            endif
            if(nxy_factor * nxy(igtyp) .gt. maxc) then
               write(*,*)'Packed output array: ',nxy(igtyp)
               write(*,*)'Exceeds dimensions: ',maxc,' ,reallocating'
               deallocate(cvar,stat=ierr)
               maxc=(nxp(igtyp)+1)*(nyp(igtyp)+1) + 1000
               maxc = nxy_factor * maxc
               allocate(cvar(maxc),stat=ierr)
               if (ierr .ne. 0) stop 'Reallocating cvar...'
            end if

            if (combine_grids) then
               FNAME_grid='MASS'
            else
               if (igtyp .eq. p_mass) FNAME_grid='MASS'
               if (igtyp .eq. p_stax) FNAME_grid='STAX'
               if (igtyp .eq. p_stay) FNAME_grid='STAY'
            endif

            if (.not. have_init(igtyp_out)) then
               fname = 'CFG_' // fname_grid // '.WRF'
               if (iprint .ge. 2) write (*,'(a,i10,1x,a)') &
                    & 'Calling pakset with unit, file=',iudata,trim(fname)
               call pakset_16(iudata,FNAME,1,nxp_pakset,nyp_pakset,nzp_pakset,use_16)
               if (combine_grids) then
                  if (nxp_pakset .ne. nxp(igtyp_out)+1 .or. nyp_pakset .ne. nyp(igtyp_out)+1) &
                     write (*,*) 'Inconsistent nxp/nyp from pakset:',&
                          nxp_pakset, nxp(igtyp_out)+1, nyp_pakset, nyp(igtyp_out)+1
               else
                  if (nxp_pakset .ne. nxp(igtyp_out) .or. nyp_pakset .ne. nyp(igtyp_out)) &
                     write (*,*) 'Inconsistent nxp/nyp from pakset:',&
                          nxp_pakset, nxp(igtyp_out), nyp_pakset, nyp(igtyp_out)
               endif
               !        standard file name for output
               fname = 'DATA_' // fname_grid // '.WRF'

!        open output data set and initialize to missing
               if (iprint .ge. 2) write (*,'(a,i10,1x,a,4i10)') &
                    & 'Calling datini with unit, file, iyr,imo,ida,ihr=', iudata,trim(fname), &
                    & IYR,IMO,IDA,IHR
               nxy_arg=nxy_factor * NXY(igtyp)
               call datini(iudata,fname,CVAR,nxy_arg,NLVL(igtyp),NVAR(1,igtyp),IG(igtyp), &
                    & IYR,IMO,IDA,IHR)
               if (iprint .ge. 1) write(*,'(2a)') 'Initialized output data set: ',trim(fname)
            endif
            have_init(igtyp) = .TRUE.
         endif !have_init(igtyp)
         
         ILTYP_LOOP: do iltyp=1,nlevtyp

            nlev_in = -99
            if (iltyp .eq. p_sfc) iladd_save = 0
            if (iltyp .eq. p_full) iladd_save = 0
            if (iltyp .eq. p_half) iladd_save = 1

            IVAR_LOOP: do ivar=1,nvars(iltyp,igtyp)
               !     read variables into 4d array (3rd and/or 4th dimension may be 1)
               !     perform any required units conversion
               
               iladd = iladd_save
               call setvar(iprint .ge. 8,handle,vlabels(ivar,iltyp,igtyp),varid,n3d,ndim,dimlen,ierr)
               if (ierr .ne. 0) cycle IVAR_LOOP
               nlev_var = 1
               if (ndim .eq. 4) nlev_var=dimlen(3)
               if (dimlen(1) .ne. nxp(igtyp) .or. dimlen(2) .ne. nyp(igtyp) .or. &
                    & (nlev_var .ne. nsigl(iltyp)-iladd .and. nlev_var .ne. nsoil)) then
                  write (*,*) 'Bad dimensions, skipping ',trim(vlabels(ivar,iltyp,igtyp))
                  cycle IVAR_LOOP
               endif
               if (nlev_var .ne. nsigl(iltyp)-iladd) iladd = -iladd_soil
               
               if (nlev_in .ne. nlev_var) then
                  if (iprint .ge. 8) write (*,*) 'Allocating arrays, previous nlev_in= ',nlev_in
                  nlev_in = nlev_var
                  if (   (iltyp .eq. p_sfc .and. &
                       &  nlev_in .ne. nsoil .and. nlev_in .ne. 1) .or. &
                       & (iltyp .eq. p_full .and. &
                       &  nlev_in .ne. nsigl(p_full)) .or. &
                       & (iltyp .eq. p_half .and. &
                       &  nlev_in .ne. nsigl(p_full)-1) ) then
                     write (*,*) 'Bad levels: iltyp, nlev_in, nsoil, nsigl: ',iltyp, nlev_in, nsoil, nsigl
                     write (*,*) '   ==> Skipping : ',trim(vlabels(ivar,iltyp,igtyp))
                     nlev_in = -99
                     cycle IVAR_LOOP
                  endif
                     
                  if (igtyp .eq. p_mass) then
                     deallocate(rvar_mass,stat=ierr)
                     allocate(rvar_mass(nxp(igtyp),nyp(igtyp),nlev_in,ntimes),stat=ierr)
                     if (ierr .ne. 0) stop 'Allocate rvar_mass...'
                  elseif (igtyp .eq. p_stax) then
                     deallocate(rvar_stax,stat=ierr)
                     allocate(rvar_stax(nxp(igtyp),nyp(igtyp),nlev_in,ntimes),stat=ierr)
                     if (ierr .ne. 0) stop 'Allocate rvar_stax...'
                  elseif (igtyp .eq. p_stay) then
                     deallocate(rvar_stay,stat=ierr)
                     allocate(rvar_stay(nxp(igtyp),nyp(igtyp),nlev_in,ntimes),stat=ierr)
                     if (ierr .ne. 0) stop 'Allocate rvar_stay...'
                  endif
               endif

               if (igtyp .eq. p_mass) then
                  status = nf_get_var_real(handle, varid, rvar_mass)
                  if(status .eq. nf_noerr .and. &
                       & (vcnvrt(ivar,iltyp,igtyp) .ne. 1.0 .or. vadd(ivar,iltyp,igtyp) .ne. 0.0)) &
                       & call datcnv(rvar_mass,nxp(igtyp),nyp(igtyp),nlev_var,ntimes,itime_local, &
                       & vcnvrt(ivar,iltyp,igtyp),vadd(ivar,iltyp,igtyp))
               elseif (igtyp .eq. p_stax) then
                  status = nf_get_var_real(handle, varid, rvar_stax)
                  if(status .eq. nf_noerr .and.  &
                       & (vcnvrt(ivar,iltyp,igtyp) .ne. 1.0 .or. vadd(ivar,iltyp,igtyp) .ne. 0.0)) &
                       & call datcnv(rvar_stax,nxp(igtyp),nyp(igtyp),nlev_var,ntimes,itime_local, &
                       & vcnvrt(ivar,iltyp,igtyp),vadd(ivar,iltyp,igtyp))
               elseif (igtyp .eq. p_stay) then
                  status = nf_get_var_real(handle, varid, rvar_stay)
                  if(status .eq. nf_noerr .and.  &
                       & (vcnvrt(ivar,iltyp,igtyp) .ne. 1.0 .or. vadd(ivar,iltyp,igtyp) .ne. 0.0)) &
                       & call datcnv(rvar_stay,nxp(igtyp),nyp(igtyp),nlev_var,ntimes,itime_local, &
                       & vcnvrt(ivar,iltyp,igtyp),vadd(ivar,iltyp,igtyp))
               endif
               if (status .ne. nf_noerr) then
                  write (*,*) 'Error getting variable=',trim(vlabels(ivar,iltyp,igtyp))
                  cycle IVAR_LOOP
               end if
               
               ILEV_LOOP: do ilev=1+iladd_save,nsigl(iltyp)
                  krec(igtyp) = krec(igtyp) + 1
                  krec_out(igtyp_out) = krec_out(igtyp_out) + 1
                  level=sigl(ilev,iltyp)
                  ! determine ARL level number
                  if (iltyp .eq. p_sfc) then
                     kl=1
                  else
                     if (.not. combine_levs .or. iltyp .eq. p_full) then
                        do kl=1,nlvl(igtyp_out)
                           if (level .eq. lvl(kl,igtyp_out)) exit
                        enddo
                     else
                        kl = ilev !half-levels are assigned to next higher full level if combine_levs
                     endif
                  endif
                  if (iprint .ge. 9) write (*,'(10(a,i10))') ' level=',sigl(ilev,iltyp), &
                       & ' igtyp=',igtyp,' iltyp=',iltyp,' ivar=',ivar,' ilev=',ilev,' kl=',kl

                  !     [if combining grids, copy into output array]
                  !     then pack into ARL format and continue
                  if (iprint .ge. 8) write (*,'(a,i10,8i10)') &
                       & 'Calling pakrec with unit, iyr,imo,ida,ihr,imn,ifh,kl,0=', iudata, &
                       & IYR,IMO,IDA,IHR,IMN,IFH,KL,0
                  if (igtyp .eq. p_mass) then
                     if (.not. combine_grids) then
                        nxy_arg=nxy_factor * NXY(igtyp)
                        call pakrec_16(iudata,rvar_mass(:,:,ilev-iladd,itime_local), &
                             & cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                             & vnames(ivar,iltyp,igtyp), &
                             & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
                     else
                        call copy_arr(rvar_mass(:,:,ilev-iladd,itime_local),nxp(igtyp),nyp(igtyp), &
                             & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
                     endif
                  elseif (igtyp .eq. p_stax) then
                     if (.not. combine_grids) then
                        nxy_arg=nxy_factor * NXY(igtyp)
                        call pakrec_16(iudata,rvar_stax(:,:,ilev-iladd,itime_local),cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                             & vnames(ivar,iltyp,igtyp), &
                             & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
                     else
                        call copy_arr(rvar_stax(:,:,ilev-iladd,itime_local),nxp(igtyp),nyp(igtyp), &
                             & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
                     endif
                  elseif (igtyp .eq. p_stay) then
                     if (.not. combine_grids) then
                        nxy_arg=nxy_factor * NXY(igtyp)
                        call pakrec_16(iudata,rvar_stay(:,:,ilev-iladd,itime_local),cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                             & vnames(ivar,iltyp,igtyp), &
                             & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
                     else
                        call copy_arr(rvar_stay(:,:,ilev-iladd,itime_local),nxp(igtyp),nyp(igtyp), &
                             & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
                     endif
                  endif
                  if (combine_grids) then
                     nxy_arg=nxy_factor * NXY(igtyp_out)
                     call pakrec_16(iudata,rvar_out,cvar,nxp(igtyp_out)+1,nyp(igtyp_out)+1,nxy_arg, &
                          & vnames(ivar,iltyp,igtyp), &
                          & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
                  end if
               end do ILEV_LOOP
            end do IVAR_LOOP
         end do ILTYP_LOOP
      end do IGTYP_LOOP
      
      do igtyp=1,ngridtyp
         write (*,'(a,i10,a,2i10)') 'For grid type ',igtyp, &
              & ' the number of processed records is:', krec(igtyp)
         write (*,'(a,i10,a,2i10)') 'For grid type ',igtyp, &
              & ' the number of processed output records is:', krec_out(igtyp)
         if (nlvl(igtyp) .gt. 0) then
            if (sum(nvar(1:nlvl(igtyp),igtyp)) .ne. krec_out(igtyp)) &
                 & write (*,'(a,i10)') &
                 & 'WARNING: The expected number of output records is:', &
                 & sum(nvar(1:nlvl(igtyp),igtyp)), &
                 & '      in the case of missing sfc variables:', &
                 & sum(nvar(2:nlvl(igtyp),igtyp))
         else
            if (krec_out(igtyp) .gt. 0) &
                 & write (*,'(a,i10)') 'WARNING: The expected number of output records is:', 0
         end if
      end do
      RETURN
      end SUBROUTINE xtract
    end module xtract_module
    
