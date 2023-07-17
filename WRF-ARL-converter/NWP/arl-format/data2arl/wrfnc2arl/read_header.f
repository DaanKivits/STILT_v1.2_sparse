! Use emacs f90-mode:   -*- mode:f90  -*-
module read_header_module

! $Id: read_header.f,v 1.1 2007/02/07 19:46:37 trn Exp $

contains    
      subroutine read_header(iu,fname,ndatlin,ierr,iprint)

      implicit none

      integer, intent(in) :: iu, iprint  !unit number, printout ctl
      character (len=*) :: fname !file name
      integer, intent(out) :: ndatlin, ierr !no of datalines, error code

      integer :: status, iline
      character(len=256) :: line

      ndatlin=0
      ierr = -1
      if (iprint .ge. 1) write (*,'(2a)') 'Reading header from file: ',trim(fname)
      open(unit=iu,file=fname,iostat=status)
      if (status .ne. 0) goto 900

      ierr = 1
      read(iu,*,iostat=status) ndatlin
      if (status .ne. 0) goto 900
      if (iprint .ge. 4) write (*,'(a,i10)') 'Decoded ndatlin= ',ndatlin

      LINE_LOOP: do iline=2,10000
         ierr = iline
         read (iu,'(a)',iostat=status) line
         if (status .ne. 0) goto 900
         if (iprint .ge. 4) write (*,'(2a)') 'Header line: ',trim(line)
         if (line(1:9) .eq. 'ENDHEADER') then
            ierr = 0
            exit LINE_LOOP
         end if
      end do LINE_LOOP
      
900   if (ierr .ne. 0 .and. iprint .ge. 1) then
         write (*,'(a,i10)') 'Error reading header with error code: ', &
              & ierr
         write (*,'(a,i10)') '(-1: opening file; >0 : line number of error)'
      end if
      return
      end subroutine read_header
    end module read_header_module
    
