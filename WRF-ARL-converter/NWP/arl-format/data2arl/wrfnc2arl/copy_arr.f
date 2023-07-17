! Use emacs f90-mode:   -*- mode:f90  -*-
module copy_arr_module

! $Id: copy_arr.f,v 1.1 2007/02/07 19:46:36 trn Exp $

contains    
      subroutine copy_arr(xin,nin1,nin2,xout,nout1,nout2)

        implicit none 
        integer, intent(in) :: nin1,nin2,nout1,nout2
        real, intent(in) :: xin(nin1,nin2)
        real, intent(out) :: xout(nout1,nout2)

        integer :: i,j

        do j=1,min(nin2,nout2)
           do i=1,min(nin1,nout1)
              xout(i,j) = xin(i,j)
           end do
        end do
! Fill extra column(s), if any
        do j=nin2+1,nout2
           do i=1,min(nin1,nout1)
              xout(i,j) = xin(i,nin2)
           end do
        end do
! Fill extra row(s), if any
        do j=1,min(nin2,nout2)
           do i=nin1+1,nout1
              xout(i,j) = xin(nin1,j)
           end do
        end do
! Fill upper right corner, if needed
        do j=nin2+1,nout2
           do i=nin1+1,nout1
              xout(i,j) = xin(nin1,nin2)
           end do
        end do
        
        return
      end subroutine copy_arr
    end module copy_arr_module
    
