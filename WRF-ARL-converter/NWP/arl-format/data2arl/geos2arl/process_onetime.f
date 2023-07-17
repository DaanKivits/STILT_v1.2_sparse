! Use emacs f90-mode:   -*- mode:f90  -*-
module process_onetime_module
  implicit none
  private
  public :: process_onetime
contains
  subroutine process_onetime(iudata, iucfg, iprint, &
       nlevtyp, p_sfc, p_lay, top_sigma, &
       nvars, ndim1v, &
       vlabels, lchar, vnames, vcnvrt, vadd, vfnames, vi_tavg, &
       mdlid, use_16, geos_prefix, geos_sysver, valid_time, geos_suffix)

    use get_fname_module
    use make_index_module
    use set_lev_module
    use set_map_module

    use datcnv_module !wrfnc2arl module
    use datini_module !wrfnc2arl module
    use getdim_module !wrfnc2arl module
    use setvar_module !wrfnc2arl module

    implicit none

    include 'netcdf.inc'
    include 'pakset_16.INC'

    integer, intent(in) :: iudata, iucfg, iprint, nlevtyp, p_sfc, p_lay, ndim1v, lchar
    integer, intent(in) :: nvars(nlevtyp)
    character (len=lchar), intent(in) :: vlabels(ndim1v,nlevtyp), vfnames(ndim1v,nlevtyp)
    character (len=4), intent(in) :: vnames(ndim1v,nlevtyp)
    real, intent(in) :: vcnvrt(ndim1v,nlevtyp),vadd(ndim1v,nlevtyp)
    integer, intent(in) :: vi_tavg(ndim1v,nlevtyp)
    character (len=4), intent(in) :: mdlid
    logical, intent(in) :: use_16
    character (len=lchar), intent(in) :: geos_prefix, geos_sysver, valid_time, geos_suffix
    real, intent(in) :: top_sigma

    integer :: handle, handle2
    character (len=lchar) :: handle_fname, handle2_fname, fname
    integer :: yr_valid, mo_valid, dy_valid, hr_valid, mi_valid, sc_valid, itime, ntimes, yy_valid
    logical :: pass_valid
    integer :: ifh, ierr
    integer :: ivar, ivar_apre, ivar_aprs, delta_per_loop, abs_delta_per_loop, minutes_unit_factor
    integer :: nx, ny, nlay, ndim1l, nxy_factor, nxy_arg, maxc, nlay_keep
    real :: grids(12), dummy(3)
    integer, allocatable :: sigl(:,:), nvar(:), lvl(:)
    integer :: nsigl(nlevtyp), iloop, nloop, nlvl, ilev, nlev_in, nlev_var
    logical :: inst_to_avg
    real    :: loop_weight, sum_weight
    integer :: ilev_bottom_up, ilev_top_down
    integer :: ndim,dimlen(nf_max_var_dims),n3d,varid
    integer :: ig, krec, krec_out, level, kl, nzp, iltyp, iladd, iladd_save

    !     packed output array
    CHARACTER (len=1), allocatable ::  CVAR(:)
    ! data arrays
    REAL, allocatable :: rvar_4d(:,:,:,:), rvar_tmp(:,:,:,:)

    krec = 0
    krec_out = 0

    ! Input: Time information: from valid_time
    call get_time(valid_time, yr_valid, mo_valid, dy_valid, hr_valid, mi_valid, sc_valid, pass_valid)
    yy_valid = mod(yr_valid,100) ! 2-digit year
    ifh=0
    if (.not. pass_valid) then
       write (*,*) 'process_onetime: bad valid_time: ',valid_time
       stop 'process_onetime: bad valid_time'
    end if
    ! Set delta_time per loop:
    if (len_trim(valid_time) .eq. 15) then
       !delta_time in seconds
       minutes_unit_factor = 60
    elseif (len_trim(valid_time) .eq. 13) then
       !delta_time in minutes
       minutes_unit_factor = 1
    else
       write (*,*) 'process_onetime: bad valid_time len_trim: ',valid_time,len_trim(valid_time)
       stop 'process_onetime: bad valid_time len_trim'
    end if
    abs_delta_per_loop=60*minutes_unit_factor ! average one-hourly files into 3-hour averages

    ! For 4d-arrays from netcdf file: itime=1, ntimes=1
    itime=1
    ntimes=1

    ! Obtain Level and grid information: 
    ! Find vfname containing APRE/APRS
    ivar_apre=-99
    ivar_aprs=-99
    do ivar=1,max(nvars(p_lay),nvars(p_sfc))
       if (ivar .le. nvars(p_lay)) then
          if (vnames(ivar,p_lay) .eq. 'APRE') ivar_apre=ivar
       end if
       if (ivar .le. nvars(p_sfc)) then
          if (vnames(ivar,p_sfc) .eq. 'APRS') ivar_aprs=ivar
       end if
    end do
    if (ivar_apre .eq. -99 .or. ivar_aprs .eq. -99) then
       write (*,*) 'process_onetime: Cannot find APRE/APRS in variable table, abort'
       stop 'process_onetime: Cannot find APRE/APRS in variable table, abort'
    end if
    iloop=1
    delta_per_loop=abs_delta_per_loop
    call get_fname(geos_prefix, vfnames(ivar_apre,p_lay), geos_sysver, valid_time, geos_suffix, &
         minutes_unit_factor,iloop, delta_per_loop, handle_fname, ierr)
    if (ierr .ne. 0) then
       write (*,*) 'process_onetime: Error generating file name for APRE'
       stop 'process_onetime: Error generating file name for APRE'
    end if
    ! open file
    if (iprint .ge. 4) write (*,'(2a)') 'opening APRE file: ',trim(handle_fname)
    ierr = NF_OPEN(handle_fname, nf_nowrite, handle)
    if (ierr .ne. nf_noerr) then
       write (*,*) 'Error opening ',trim(handle_fname),' : ',nf_strerror(ierr)
       stop 'process_onetime: Error opening file for APRE'
    end if
    call get_fname(geos_prefix, vfnames(ivar_aprs,p_sfc), geos_sysver, valid_time, geos_suffix, &
         minutes_unit_factor,iloop, delta_per_loop, handle2_fname, ierr)
    if (ierr .ne. 0) then
       write (*,*) 'process_onetime: Error generating file name for APRS'
       stop 'process_onetime: Error generating file name for APRS'
    end if
    if (trim(handle_fname) .eq. trim(handle2_fname)) then
       handle2 = handle
    else
       ! open file
       if (iprint .ge. 4) write (*,'(2a)') 'opening APRS file: ',trim(handle2_fname)
       ierr = NF_OPEN(handle2_fname, nf_nowrite, handle2)
       if (ierr .ne. nf_noerr) then
          write (*,*) 'Error opening ',trim(handle2_fname),' : ',nf_strerror(ierr)
          stop 'process_onetime: Error opening file for APRS'
       end if
    end if
    
    ! obtain nx(lon),ny(lat),nlay(lev) grid dimension:
    call getdim(iprint .ge. 8,handle,'lon',nx,ierr)
    if (ierr .eq. 0) call getdim(iprint .ge. 8,handle,'lat',ny,ierr)
    if (ierr .eq. 0) call getdim(iprint .ge. 8,handle,'lev',nlay,ierr)
    ! call set_map: read lon,lat, construct grids info
    if (ierr .eq. 0) call set_map(iprint .ge. 8,handle,nx,ny,grids,ierr)
    if (ierr .ne. 0) then
       write (*,*) 'process_onetime: Error obtaining dim (lon,lat,lev) or map info'
       stop 'process_onetime: Error obtaining dim (lon,lat,lev) or map info'
    end if
    if (vfnames(ivar_apre,p_lay) .ne. vfnames(ivar_aprs,p_sfc)) then
       write (*,*) 'process_onetime: APRE and APRS are not in the same file, assuming dims etc are the same'
    end if
    call set_lev(iprint,handle,vlabels(ivar_apre,p_lay),handle2,vlabels(ivar_aprs,p_sfc),&
         sigl,nsigl,nlay_keep,top_sigma,nlevtyp,p_sfc,p_lay,nx,ny,nlay,ierr)
    if (ierr .ne. 0) then
       write (*,*) 'process_onetime: set_lev error, abort'
       stop 'process_onetime: set_lev error, abort'
    end if
    ndim1l=nlay_keep+1
    allocate(nvar(2*ndim1l),lvl(2*ndim1l))

    if (iprint .ge. 2) &
         write (*,'(4(a,i10))') 'Grid dimensions: nx=',nx,', ny=',ny,', nlay=',nlay, &
         ', nlay_keep=',nlay_keep

    ierr=nf_close(handle)
    if (ierr .ne. nf_noerr) write (*,*) 'Error closing ',trim(handle_fname), &
         ' : ',nf_strerror(ierr)
    if (handle .ne. handle2) then
       ierr=nf_close(handle2)
       if (ierr .ne. nf_noerr) write (*,*) 'Error closing ',trim(handle2_fname), &
            ' : ',nf_strerror(ierr)
    end if
    
    ! allocate output arrays
    nxy_factor = 1
    if (use_16) nxy_factor = 2
    maxc=nxy_factor * (nx * ny + 1000)
    allocate(cvar(maxc))

    ! call ***make_index*** : write cfg info to ascii file

    fname = 'CFG_' // trim(valid_time) // '.GEOS'
    if (iprint .ge. 2) write (*,'(a,i10,1x,a)') &
         & 'Writing cfg info to unit, file=',iucfg,trim(fname)
    OPEN(iucfg,file=trim(fname))

    call make_index(nlevtyp, p_sfc, p_lay, &
         vnames, nvars, nsigl, ndim1l, ndim1v, &
         sigl, iprint, iucfg, nx, ny, mdlid, grids, &
         nvar, lvl, nlvl)

    CLOSE (iucfg)

    ! call ***pakset_16***: read cfg info, initialize data structures
    if (iprint .ge. 8) write (*,'(a,i10,1x,a)') &
         & 'Calling pakset with unit, file=',iudata,trim(fname)
    call pakset_16(iudata,FNAME,1,nx,ny,nzp,use_16)
    !        standard file name for output
    fname = 'DATA_' // trim(valid_time) // '.GEOS'

    ! call ***datini***: open output data set and initialize to missing
    if (iprint .ge. 3) write (*,'(a,i10,1x,a,4i10)') &
         & 'Calling datini with unit, file, iyr,imo,ida,ihr=', iudata,trim(fname), &
         & yy_valid, mo_valid, dy_valid, hr_valid
    nxy_arg=nxy_factor * NX * ny
    ig=99
    call datini(iudata,fname,CVAR,nxy_arg,NLVL,NVAR(1),ig, &
         & yy_valid, mo_valid, dy_valid, hr_valid)
    if (iprint .ge. 1) write(*,'(2a)') 'Initialized output data set: ',trim(fname)

    ! Loop over level type:
    ! Initialize: nlev_in, set iladd_save (1 for layer, 0 otherwise)
    ILTYP_LOOP: do iltyp=1,nlevtyp

       nlev_in = -99
       if (iltyp .eq. p_sfc) iladd_save = 0
       if (iltyp .eq. p_lay) iladd_save = 1

       ! Loop over variables for current level type:
       IVAR_LOOP: do ivar=1,nvars(iltyp)
          !     read variables into 4d array (3rd and/or 4th dimension may be 1)
          !     perform any required units conversion

          !  reset iladd=iladd_save

          iladd = iladd_save

          !  set up time averaging as needed:

          if (vi_tavg(ivar,iltyp) .ge. 0) then
             !   i_tavg=0 (instantaneous valid at end of time period), 1 (average over output interval): 
             !            read the current file, delta t=0
             !   i_tavg>1 (average over a shorter interval): 
             !            read and average for delta t = 0, -1, ...
             nloop=max(1,vi_tavg(ivar,iltyp))
             delta_per_loop=-abs_delta_per_loop
             inst_to_avg = .FALSE.
          elseif (vi_tavg(ivar,iltyp) .eq. -1) then
             !   i_tavg=-1 (average over a shorter interval for instantaneous output values): 
             !             read and average for delta t = 0, +1
             nloop=2
             delta_per_loop=abs_delta_per_loop
             inst_to_avg = .FALSE.
          elseif (vi_tavg(ivar,iltyp) .lt. -1) then
             !   i_tavg<-1 compute average from instantaneous values
             !             read and average for delta t = 0, ...,  -navg (use 1/2 weight at endpoints)
             nloop=-vi_tavg(ivar,iltyp)
             inst_to_avg = .TRUE.
             ! inst3 files are spaced apart the same distance as tavg3 files, need to use 3-hour file spacing
             delta_per_loop=-3*abs_delta_per_loop
          else
             write (*,'(a,i10,4a)') 'Invalid vi_tavg=',vi_tavg(ivar,iltyp),', skipping vname=',vnames(ivar,iltyp), &
                  ' ,vlabel=',trim(vlabels(ivar,iltyp))
             cycle IVAR_LOOP
          end if

          sum_weight=0.
          ITIM_LOOP: do iloop=1,nloop
             ! For each file (time period): set dt, weight, call ***get_fname***
             ! Note: in get_fname, file time stamp is determined from valid_time+(iloop-1)*dt
             if (inst_to_avg .and. (iloop .eq. 1 .or. iloop .eq. nloop)) then
                loop_weight=0.5 ! use half-weight at endpoint
             else
                loop_weight=1.0
             end if
             
             call get_fname(geos_prefix, vfnames(ivar,iltyp), geos_sysver, valid_time, geos_suffix, &
                  minutes_unit_factor, iloop, delta_per_loop, handle_fname, ierr)
             if (ierr .ne. 0) then
                write (*,'(4a)') 'process_onetime: Error generating file name, skipping vname=',vnames(ivar,iltyp), &
                     ' ,vlabel=',trim(vlabels(ivar,iltyp))
                cycle IVAR_LOOP
             end if
             if (iprint .ge. 3) &
                  write (*,'(a,2i5,3a,i5,a,f6.2,2a)') 'Processing iltyp,ivar=',iltyp,ivar,': ',trim(vlabels(ivar,iltyp)), &
                  ', iloop=',iloop,', weight=',loop_weight,', fname=',trim(handle_fname)
             ! open file, call ***setvar*** for vlabels(ivar,iltyp)
             !  skip to end of var_loop if file or variable missing
             ierr = NF_OPEN(handle_fname, nf_nowrite, handle)
             if (ierr .ne. nf_noerr) then
                write (*,'(3a)') 'Error opening ',trim(handle_fname),' : ',nf_strerror(ierr)
                cycle IVAR_LOOP
             end if
             call setvar(iprint .ge. 8,handle,vlabels(ivar,iltyp),varid,n3d,ndim,dimlen,ierr)
             if (ierr .ne. 0) then
                write (*,'(a,2i10,3a,i10,2a)') 'Skipping variable because of setvar error for iltyp,ivar=', &
                     iltyp,ivar,': ',trim(vlabels(ivar,iltyp)), &
                     ', iloop=',iloop,', fname=',trim(handle_fname)
                cycle IVAR_LOOP
             end if
             !    set nlev_var, check dimensions against nxp, nyp, 
             !    nsigl-iladd for iltyp: skip if mismatch
             nlev_var = 1
             if (ndim .eq. 4) nlev_var=dimlen(3)
             if (dimlen(1) .ne. nx .or. dimlen(2) .ne. ny .or. &
                  & (nlev_var .lt. nsigl(iltyp)-iladd)) then
                write (*,'(2a)') 'Bad dimensions, skipping ',trim(vlabels(ivar,iltyp))
                cycle IVAR_LOOP
             endif

             !  Allocate if needed (nlev_in .ne. nlev_var):
             !   4d array rvar_4d: nx,ny,nlev_in,ntimes
             if (nlev_in .ne. nlev_var) then
                if (iprint .ge. 8) write (*,'(2(a,i10))') 'Allocating arrays with nlev_var=',nlev_var,&
                     ', previous nlev_in= ',nlev_in
                nlev_in = nlev_var
                if ( (iltyp .eq. p_sfc .and. nlev_in .ne. 1) .or. &
                     (iltyp .eq. p_lay .and. nlev_in .lt. nsigl(p_lay)-1) ) then
                   write (*,'(a,3i10)') 'Bad levels: iltyp, nlev_in, nsigl: ',iltyp, nlev_in, nsigl
                   write (*,'(2a)') '   ==> Skipping : ',trim(vlabels(ivar,iltyp))
                   nlev_in = -99
                   cycle IVAR_LOOP
                endif
                deallocate(rvar_4d,rvar_tmp,stat=ierr)
                allocate(rvar_4d(nx,ny,nlev_in,ntimes),rvar_tmp(nx,ny,nlev_in,ntimes),stat=ierr)
                if (ierr .ne. 0) stop 'Allocate rvar_4d...'
             end if
             !    call (***nf_get_var_real***), 
             ierr = nf_get_var_real(handle, varid, rvar_tmp)
             ! close file
             ierr = NF_close(handle)
             if (ierr .ne. nf_noerr) then
                write (*,'(3a)') 'Error closing ',trim(handle_fname),' : ',nf_strerror(ierr)
                cycle IVAR_LOOP
             end if
             !   apply averaging
             sum_weight=sum_weight+loop_weight
             if (iloop .eq. 1) then
                rvar_4d=rvar_tmp*loop_weight
             else
                rvar_4d=rvar_4d + rvar_tmp*loop_weight
             end if
             if (iloop .eq. nloop) then
                rvar_4d=rvar_4d / sum_weight
                !    special handling for delp: pass in dummy aprs, set diag to false:
                if (vlabels(ivar,iltyp) .eq. 'DELP') &
                     call delp2pre(.FALSE.,rvar_4d,rvar_tmp,nx,ny,nlay,dummy)
                ! call ***datcnv**: rescale variable
                if (vcnvrt(ivar,iltyp) .ne. 1.0 .or. vadd(ivar,iltyp) .ne. 0.0) &
                     call datcnv(rvar_4d,nx,ny,nlev_var,ntimes,itime, &
                     vcnvrt(ivar,iltyp),vadd(ivar,iltyp))
             end if
          end do ITIM_LOOP

          !  Loop over levels: ilev=1+iladd_save, nsigl(iltyp)
          ILEV_LOOP: do ilev=1+iladd_save,nsigl(iltyp)
             !   Increment output record counter
             krec = krec + 1
             krec_out = krec_out + 1
             !   Set level value: level=sigl(ilev) and 
             level=sigl(ilev,iltyp)
             !index into rvar for bottom_up storage of rvar: 1 for sfc, ilev-1 for lay
             ilev_bottom_up = ilev-iladd
             ! determine ARL level number and index into rvar
             if (iltyp .eq. p_sfc) then
                ! for surface (2d) variables:
                !  ilev=1 (iladd_save=0), iladd=0, ilev-iladd=1, kl=1
                kl=1
                !index into rvar for sfc: 1
                ilev_top_down = 1
             else
                ! for layer values:
                !  ilev=2:nsigl(p_lay) (iladd_save=1), ilev-iladd=ilev-1=1:no of half levels
                !  kl=ilev=2:nsigl(p_lay) - recall that nsigl(p_lay)=1+no of half levels
                kl = ilev
                !index into rvar for top_down storage of rvar: (nlay,nlay-nlay_keep+1,-1)
                ilev_top_down = nlay - ilev_bottom_up + 1
             endif

             if (iprint .ge. 8) write (*,'(10(a,i10))') ' level=',sigl(ilev,iltyp), &
                   ' iltyp=',iltyp,' ivar=',ivar,' ilev=',ilev,' kl=',kl, &
                   ' ilev_top_down=',ilev_top_down

             ! then pack into ARL format and continue: call ***pakrec_16***, writes data output record

             if (iprint .ge. 8) write (*,'(a,i10,1x,a,8i10)') &
                  & 'Calling pakrec with unit, iyr,imo,ida,ihr,imn,ifh,kl,0=', iudata, vnames(ivar,iltyp), &
                  & yy_valid,mo_valid,dy_valid,hr_valid,mi_valid,IFH,KL,0
             call pakrec_16(iudata,rvar_4d(:,:,ilev_top_down,itime), &
                  cvar,nx,ny,nxy_arg, &
                  vnames(ivar,iltyp), &
                  yy_valid,mo_valid,dy_valid,hr_valid,mi_valid,IFH,KL,0,use_16)
          end do ILEV_LOOP
       end do IVAR_LOOP
    end do ILTYP_LOOP

    if (iprint .ge. 1) write (*,'(a,i10,a,2i10)') 'The number of processed records is:', krec
    if (iprint .ge. 1) write (*,'(a,i10,a,2i10)') 'The number of processed output records is:', krec_out
    if (nlvl .gt. 0) then
       if (sum(nvar(1:nlvl)) .ne. krec_out) &
            write (*,'(a,i10)') &
            'WARNING: The expected number of output records is:', &
            sum(nvar(1:nlvl)), &
            '      in the case of missing sfc variables:', &
            sum(nvar(2:nlvl))
    else
       if (krec_out .gt. 0) &
            & write (*,'(a,i10)') 'WARNING: The expected number of output records is:', 0
    end if

  end subroutine process_onetime
end module process_onetime_module
