! Use emacs f90-mode:   -*- mode:f90  -*-
module set_lev_module

  implicit none
  private
  public :: set_lev, delp2pre

contains

  subroutine set_lev(iprint,handle,vlabel_apre,handle2,vlabel_aprs,sigl,nsigl,nlay_keep,&
       top_sigma,nlevtyp,p_sfc,p_lay,nx,ny,nlay,ierr)

    use setvar_module

    implicit none
    include 'netcdf.inc'

    integer, intent(in) :: iprint
    integer, intent(in) :: handle, handle2, nlevtyp, p_sfc, p_lay, nx, ny, nlay
    character (len=*), intent(in) :: vlabel_apre, vlabel_aprs
    integer, intent(out) :: nsigl(nlevtyp), ierr, nlay_keep
    integer, allocatable, intent(out) :: sigl(:,:)
    real, intent(in) :: top_sigma

    real :: apre(nx,ny,nlay,1), aprs(nx,ny,1), sigma(nx,ny), sigma_mean, si1stats(3)
    integer :: varid, n3d, ndim, dimlen(nf_max_var_dims),k,k_bottom_up,sigma_lay(nlay+1)

    !  read in pre
    call setvar(iprint .ge. 8,handle,vlabel_apre,varid,n3d,ndim,dimlen,ierr)
    if (ierr .eq. 0) then
       if (ndim .ne. 4 .or. nx .ne. dimlen(1) .or. ny .ne. dimlen(2) .or. &
            nlay .ne. dimlen(3) .or. 1 .ne. dimlen(4)) ierr=1
    end if
    if (ierr .ne. 0) then
       write (*,'(2a)') 'set_lev: setvar error for apre: ',trim(vlabel_apre)
       return
    end if
    ierr = nf_get_var_real(handle,varid,apre)
    if (ierr .ne. nf_noerr) then
       write (*,*) 'set_lev: error getting apre values'
       return
    end if
    
    !  read in prs
    call setvar(iprint .ge. 8,handle2,vlabel_aprs,varid,n3d,ndim,dimlen,ierr)
    if (ierr .eq. 0) then
       if (ndim .ne. 3 .or. nx .ne. dimlen(1) .or. ny .ne. dimlen(2) .or. &
            1 .ne. dimlen(3)) ierr=1
    end if
    if (ierr .ne. 0) then
       write (*,'(2a)') 'set_lev: setvar error for aprs: ',trim(vlabel_aprs)
       return
    end if
    ierr = nf_get_var_real(handle2,varid,aprs)
    if (ierr .ne. nf_noerr) then
       write (*,*) 'set_lev: error getting aprs values'
       return
    end if
    
    !  special handling: if label='delp', convert to layer pressures (check against aprs)
    si1stats(:)=1.0
    if (vlabel_apre .eq. 'DELP') &
         call delp2pre(iprint .ge. 3,apre,aprs,nx,ny,nlay,si1stats)
    
!  compute mean si, store in sigl(2:nsigl,p_half)
! Note: 
!  nsigl(p_sfc)=1; sigl(1,p_sfc)=10000
!  nsigl(p_lay)=1+no of layers; sigl(1,p_lay)=10000; sigl(2:nsigl,p_lay)=si
    sigma_lay(1) = 10000
    if (iprint .ge. 3) write (*,'(2a8,a10,3a15,/,2i8,i10,3e15.6)') &
         'k_lay','k_sig','sigma','mean p/ps','min p/ps','max p/ps', &
         0,1,sigma_lay(1),si1stats

    ! Loop over layers, bottom to top (data is stored top-down)
    nlay_keep=-99
    do k=nlay,1,-1
       k_bottom_up = nlay-k+1 ! index for bottom-up storage of output
       sigma(:,:)=apre(:,:,k,1)/aprs(:,:,1)
       sigma_mean=sum(sigma)/(nx*ny)
       sigma_lay(k_bottom_up+1)=nint(10000.*sigma_mean)
       if (sigma_mean .ge. top_sigma) then
          if (iprint .ge. 3) write (*,'(2i8,i10,3e15.6)') &
               k,k_bottom_up+1,sigma_lay(k_bottom_up+1),sigma_mean,minval(sigma),maxval(sigma)
       else
          if (iprint .ge. 3) write (*,'(i8,a8,i10,3e15.6)') &
               k,' - ',sigma_lay(k_bottom_up+1),sigma_mean,minval(sigma),maxval(sigma)
          if (nlay_keep .le. 0) nlay_keep=nlay-k
       end if
    end do
    if (nlay_keep .le. 0) nlay_keep=nlay
    nsigl(p_sfc)=1
    nsigl(p_lay)=nlay_keep+1
    allocate(sigl(nlay_keep+1,nlevtyp))
    sigl(1,p_sfc)=10000
    sigl(1,p_lay)=10000
    sigl(2:(nlay_keep+1),p_lay) = sigma_lay(2:(nlay_keep+1))
    return
  end subroutine set_lev

  subroutine delp2pre(diag,apre,aprs,nx,ny,nlay,si1stats)
    implicit none
    logical, intent(in) :: diag
    integer, intent(in) :: nx, ny, nlay
    real, intent(inout) :: apre(nx,ny,nlay,1)
    real, intent(in) :: aprs(nx,ny,1)
    real, intent(out) :: si1stats(3)

    real :: ptop=1. ! 0.01 mb
    real :: sigma(nx,ny)
    integer :: k

    if (diag) write (*,'(a,g15.6)') 'Converting delp to p using ptop=',ptop
    ! First convert to edge pressures:
    apre(:,:,1,1) = ptop + apre(:,:,1,1)
    do k=2,nlay
       apre(:,:,k,1) = apre(:,:,k-1,1) + apre(:,:,k,1)
    end do
    
    ! Diagnostic: check bottom against aprs
    if (diag) then
       sigma(:,:) = apre(:,:,nlay,1)/aprs(:,:,1)
       si1stats(1) = sum(sigma)/(nx*ny)
       si1stats(2) = minval(sigma)
       si1stats(3) = maxval(sigma)
    endif
    ! Now convert to layer pressure
    do k=nlay,2,-1
       apre(:,:,k,1) = 0.5*(apre(:,:,k-1,1)+apre(:,:,k,1))
    end do
    apre(:,:,1,1) = 0.5*apre(:,:,1,1) + 0.5*ptop
    return
  end subroutine delp2pre

end module set_lev_module
