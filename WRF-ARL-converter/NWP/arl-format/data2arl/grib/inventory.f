      PROGRAM INVENTORY

!-------------------------------------------------------------------------------
! Produces and inventory listing of all the records in grib file
!-------------------------------------------------------------------------------
! Last Revised: 26 Jun 1996
!               21 Dec 2000 (RRD) - fortran90 upgrade
!               21 Dec 2001 (RRD) - ecmwf compatibility
!               18 Jun 2002 (RRD) - corrected time range indicator
!               18 Mar 2004 (RRD) - print out time range indicator
!-------------------------------------------------------------------------------

      CHARACTER fname*40, buff(150000)*1

      WRITE(*,*)'Enter file name: '
      READ(*,'(a)')fname
      OPEN(10,file=fname,form='unformatted',access='direct',recl=1)

      WRITE(*,*)'Inventory dump of file: ',fname
      WRITE(*,*)
      WRITE(*,*)'yr mo da hr  f1  f2  fx     varib   level      type grid numx numy'

      krec=1
  100 CONTINUE

!---->find grib start
      call fcread(10,buff,krec,4,kret)
      if(kret.ne.0)goto 900
      if(buff(1)//buff(2)//buff(3)//buff(4).ne.'GRIB')then
         krec=krec+1 
         goto 100
      end if

!---->indicator section (8 bytes)
      call fcread(10,buff,krec,8,kret)
      if(kret.ne.0)goto 900
      call sect0(buff,8,ktot)

!---->product definition section (unkown length - usually 28)
      krec=krec+8
      call fcread(10,buff,krec,28,kret)
      if(kret.ne.0)goto 900
      call sect1(buff,28,klen,kgrd,kvar,ktyp,level,                            &
                 iyr,imo,ida,ihr,kf1,kf2,kfx,kbms)

!---->grid description section (usually 32)
      if(kbms.eq.128.or.kbms.eq.192)then
         krec=krec+klen
         call fcread(10,buff,krec,32,kret)
         if(kret.ne.0)goto 900
         call sect2(buff,10,klen,ntyp,nipt,njpt)
      end if

!---->print out field summary
      WRITE(*,'(4i3,3i4,5x,i5,i8,5x,4i5)')                                    &
         iyr,imo,ida,ihr,kf1,kf2,kfx,kvar,level,ktyp,kgrd,nipt,njpt

!---->bit map section
      if(kbms.eq.64.or.kbms.eq.192)then
         krec=krec+klen
         call fcread(10,buff,krec,6,kret)
         if(kret.ne.0)goto 900
         call sect3(buff,6,klen)
      end if

!---->binary data section (descriptors up to 14)
      krec=krec+klen
      call fcread(10,buff,krec,14,kret)
      if(kret.ne.0)goto 900
      call sect4(buff,4,klen)
      krec=krec+klen

!---->grib 7777 end section (always 4 bytes)
      krec=krec+4
      GO TO 100

  900 CLOSE(10)
      END


      SUBROUTINE fcread(kunit,buff,krec,nbyte,kret)
      CHARACTER buff(nbyte)*1
      DO k=1,nbyte
         READ(kunit,rec=krec+k-1,iostat=kret)buff(k)
      END DO
      RETURN
      END


      SUBROUTINE SECT0(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      klen=ichar(buff(7))+ichar(buff(6))*256+ichar(buff(5))*65536
      RETURN
      END


      SUBROUTINE SECT1(buff,nbyte,klen,kgrd,kvar,ktyp,level,                   &
         iyr,imo,ida,ihr,kf1,kf2,kfx,kbms)
      CHARACTER buff(nbyte)*1
      klen=ichar(buff(3))+ichar(buff(2))*256+ichar(buff(1))*65536
      kgrd=ichar(buff(7))
      kbms=ichar(buff(8))
      kvar=ichar(buff(9))
      ktyp=ichar(buff(10))
      level=ichar(buff(12))+ichar(buff(11))*256
      iyr=ichar(buff(13))
      imo=ichar(buff(14))
      ida=ichar(buff(15))
      ihr=ichar(buff(16))
      IF(ichar(buff(21)).EQ.10)THEN
         kf2=ichar(buff(20))+ichar(buff(19))*256
         kf1=kf2
      ELSEIF(ichar(buff(18)).GE.8)THEN
         kf1=ichar(buff(19))*ichar(buff(18))
         kf2=ichar(buff(20))*ichar(buff(18))
      ELSE
         kf1=ichar(buff(19))
         kf2=ichar(buff(20))
      END IF
      kfx=ichar(buff(21))

      RETURN
      END

      SUBROUTINE SECT2(buff,nbyte,klen,ntyp,nipt,njpt)
      CHARACTER buff(nbyte)*1
      klen=ichar(buff(3))+ichar(buff(2))*256+ichar(buff(1))*65536
      ntyp=ichar(buff(6))
      nipt=ichar(buff(8))+ichar(buff(7))*256
      njpt=ichar(buff(10))+ichar(buff(9))*256
      RETURN
      END

      SUBROUTINE SECT3(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      klen=ichar(buff(3))+ichar(buff(2))*256+ichar(buff(1))*65536
      RETURN
      END

      SUBROUTINE SECT4(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      klen=ichar(buff(3))+ichar(buff(2))*256+ichar(buff(1))*65536
      RETURN
      END
