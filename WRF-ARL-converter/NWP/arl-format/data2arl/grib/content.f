      PROGRAM CONTENT

!-------------------------------------------------------------------------------
! Dumps the contents of all grib sections for each record
!-------------------------------------------------------------------------------
! Last Revised: 19 Feb 1999
!               21 Dec 2000 (RRD) - fortran90 upgrade
!               21 Dec 2001 (RRD) - ecmwf compatibility
!               13 Feb 2002 (RRD) - dump PV field of GDS
!               18 Jun 2002 (RRD) - corrected time range indicator
!               25 Feb 2003 (RRD) - ichar function replacement
!-------------------------------------------------------------------------------

      CHARACTER fname*80, buff(150000)*1

      WRITE(*,*)'Enter file name: '
      READ(*,'(a)')fname
      OPEN(10,file=fname,form='unformatted',access='direct',recl=1)

      krec=1
  100 CONTINUE

!---->find grib start
      call fcread(10,buff,krec,4,kret)
      if(kret.ne.0)goto 900
      if(buff(1)//buff(2)//buff(3)//buff(4).ne.'GRIB')then
         krec=krec+1 
         goto 100
      end if

!---->indicator "section 0" (8 bytes)
      call fcread(10,buff,krec,8,kret)
      call sect0(buff,8,klen)
      READ(*,*)

!---->product definition "section 1" (unkown length - usually 28)
      krec=krec+8
      call fcread(10,buff,krec,28,kret)
      if(kret.ne.0)goto 900
      call sect1(buff,28,klen,kbms)
      READ(*,*)

!---->grid description "section 2"
      if(kbms.eq.128.or.kbms.eq.192)then
         krec=krec+klen 
         call fcread(10,buff,krec,520,kret)
         if(kret.ne.0)goto 900
         call sect2(buff,520,klen)
         READ(*,*)
      end if

!---->bit map "section 3"
      if(kbms.eq.64.or.kbms.eq.192)then
         krec=krec+klen
         call fcread(10,buff,krec,6,kret)
         if(kret.ne.0)goto 900
         call sect3(buff,6,klen)
      else
         write(*,*)'Bit Map section length: 0'
         write(*,*)' '
      end if

!---->binary data "section 4"
      krec=krec+klen
      call fcread(10,buff,krec,11,kret)
      if(kret.ne.0)goto 900
      call sect4(buff,11,klen)
      READ(*,*)

!---->grib end "section 5"
      krec=krec+klen
      call fcread(10,buff,krec,4,kret)
      if(kret.ne.0)goto 900
      call sect5(buff,4)
      krec=krec+4
      READ(*,*)
      GO TO 100

  900 CLOSE(10)
      END

!-----------------------------------------------------------------

      SUBROUTINE fcread(kunit,buff,krec,nbyte,kret)
      CHARACTER buff(nbyte)*1
      DO k=1,nbyte
         READ(kunit,rec=krec+k-1,iostat=kret)buff(k)
      END DO
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT0(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      WRITE(*,*)'Identification:     ',buff(1),buff(2),buff(3),buff(4)
      klen=jchar(buff(7))+jchar(buff(6))*256+jchar(buff(5))*65536
      WRITE(*,*)'Length of message: ',klen
      WRITE(*,*)'Grib Edition Numb: ',jchar(buff(8))
      WRITE(*,*)
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT1(buff,nbyte,klen,kbms)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      WRITE(*,*)'Product definition section: ',klen
      WRITE(*,*)'Parameter Table Version:    ',jchar(buff(4))
      WRITE(*,*)'Center Identification:      ',jchar(buff(5))
      WRITE(*,*)'Model identification:       ',jchar(buff(6))
      WRITE(*,*)'Grid Identification:        ',jchar(buff(7))
      kbms=jchar(buff(8))
      WRITE(*,*)'GDS (128) or BMS (64):      ',kbms
      WRITE(*,*)'Parameter or Units:         ',jchar(buff(9))
      WRITE(*,*)'Type of level:              ',jchar(buff(10))
      level=jchar(buff(12))+jchar(buff(11))*256
      WRITE(*,*)'Height or pressure:         ',level

      WRITE(*,*)' '
      iyr=jchar(buff(13))
      imo=jchar(buff(14))
      ida=jchar(buff(15))
      ihr=jchar(buff(16))
      imn=jchar(buff(17))
      WRITE(*,*)'Initial time:   ',iyr,imo,ida,ihr,imn
      WRITE(*,*)'Time units:     ',jchar(buff(18))
      WRITE(*,*)'Time range:     ',jchar(buff(21))
      IF(jchar(buff(21)).EQ.10)THEN
         WRITE(*,*)'Forecast x:     ',jchar(buff(20))+jchar(buff(19))*256
      ELSE
         WRITE(*,*)'Forecast 1:     ',jchar(buff(19))
         WRITE(*,*)'Forecast 2:     ',jchar(buff(20))
      END IF
      kavg=jchar(buff(23))+jchar(buff(22))*256
      WRITE(*,*)'Average:        ',kavg
      WRITE(*,*)'Missing:        ',jchar(buff(24))
      kfact=jchar(buff(28))+jchar(buff(27))*256
!     when sign bit turned on (2^15) then negative number
      if(iand(kfact,32768).ne.0)kfact=-iand(kfact,32767)
      WRITE(*,*)'Decimal factor: ',kfact
      WRITE(*,*)
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT2(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      REAL sigma(100),offset(100)
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      WRITE(*,*)'Grid description:      ',klen
      WRITE(*,*)'Number of Levels:      ',jchar(buff(4))
      WRITE(*,*)'Octet location:        ',jchar(buff(5))
      WRITE(*,*)'Data type (table 6):   ',jchar(buff(6))
      nipt=jchar(buff(8))+jchar(buff(7))*256
      njpt=jchar(buff(10))+jchar(buff(9))*256
      WRITE(*,*)'Numb of i,j pts:       ',nipt,njpt
      WRITE(*,*)

!     ecmwf levels array (not used by ncep)
      if(jchar(buff(4)).eq.0)return
      k1=0
      k2=0
      do kv=1,jchar(buff(4))
         k=jchar(buff(5))+(kv-1)*4
         mexp=iand(jchar(buff(k)),127)
         mant=jchar(buff(k+3))+jchar(buff(k+2))*256+jchar(buff(k+1))*65536
         rval=16.0**(mexp-64)*(mant/16777216.0)
         if(iand(jchar(buff(k)),128).ne.0)rval=-rval
         if(kv.le.jchar(buff(4))/2)then
            k1=k1+1
            offset(k1)=rval/100.0
         else
            k2=k2+1
            sigma(k2)=rval
         end if
      end do
      WRITE(*,'(A)')'Level    Offset     Sigma'     
      do kk=1,k1
         WRITE(*,'(i5,f10.3,f10.6)')kk,offset(kk),sigma(kk)      
      end do
      WRITE(*,*)

!     when sign bit turned on (2^23) then negative number
      lat1=jchar(buff(13))+jchar(buff(12))*256+jchar(buff(11))*65536
      if(iand(lat1,8388608).ne.0)lat1=-iand(lat1,8388607)

      lon1=jchar(buff(16))+jchar(buff(15))*256+jchar(buff(14))*65536
      if(iand(lon1,8388608).ne.0)lon1=-iand(lon1,8388607)
      WRITE(*,*)'Grid 1,1 lat,lon: ',lat1,lon1
      WRITE(*,*)'Components:       ',jchar(buff(17))

      lorn=jchar(buff(20))+jchar(buff(19))*256+jchar(buff(18))*65536
      if(iand(lorn,8388608).ne.0)lorn=-iand(lorn,8388607)
      WRITE(*,*)'Last Point Lat:   ',lorn
      WRITE(*,*)'Long Orient:      ',lorn

      mxds=jchar(buff(23))+jchar(buff(22))*256+jchar(buff(21))*65536
      WRITE(*,*)'Grid x size (m):  ',mxds
      if(iand(mxds,8388608).ne.0)mxds=-iand(mxds,8388607)
      WRITE(*,*)'Last Point Long:  ',mxds

      myds=jchar(buff(26))+jchar(buff(25))*256+jchar(buff(24))*65536
      WRITE(*,*)'Grid y size (m):  ',myds
      WRITE(*,*)'Mercator latitude:',myds

      ldin=jchar(buff(25))+jchar(buff(24))*256
      WRITE(*,*)'Long Direc Incr:  ',ldin

      ldin=jchar(buff(27))+jchar(buff(26))*256
      WRITE(*,*)'Lat  Direc Incr:  ',ldin

      WRITE(*,*)'Projection center:',jchar(buff(27))
      WRITE(*,*)'Scanning mode:    ',jchar(buff(28))
      WRITE(*,*)
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT3(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      WRITE(*,*)'Bit Map section length:',klen
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT4(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      WRITE(*,*)'Binary data section:   ',klen
      kpack=jchar(buff(4))
      WRITE(*,*)'Table 11 flag:         ',kpack
      if(iand(kpack,128).ne.0)write(*,*)'Spherical Harmonics'
      if(iand(kpack, 64).ne.0)write(*,*)'Complex Packing'
      if(iand(kpack, 32).ne.0)write(*,*)'Original data integer'
      if(iand(kpack, 16).ne.0)write(*,*)'Octet 14 extension'

      kfact=jchar(buff(6))+jchar(buff(5))*256
      if(iand(kfact,32768).ne.0)kfact=-iand(kfact,32767)
      WRITE(*,*)'Binary scale factor:   ',kfact

!     floating point value consists of 7 bit exponent and 24 bit mantissa
      mexp=iand(jchar(buff(7)),127)
      mant=jchar(buff(10))+jchar(buff(9))*256+jchar(buff(8))*65536

!     value equals 2^(-24)*MANT*16^(MEXP-64)
      refval=16.0**(mexp-64)*mant/16777216.0
      if(iand(jchar(buff(7)),128).ne.0)refval=-refval
      WRITE(*,*)'Reference value:       ',refval

      WRITE(*,*)'Bit packing:           ',jchar(buff(11))
      WRITE(*,*)
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT5(buff,nbyte)
      CHARACTER buff(nbyte)*1
      WRITE(*,*)'End marker:    ',buff(1),buff(2),buff(3),buff(4)
      WRITE(*,*)
      RETURN
      END
