      PROGRAM UNPACKER

!-------------------------------------------------------------------------------
! Decodes all the grib records
!-------------------------------------------------------------------------------
! Last Revised: 19 Feb 1999
!               21 Dec 2000 (RRD) - fortran90 upgrade
!               21 Dec 2001 (RRD) - ecmwf compatibility
!               21 Aug 2002 (RRD) - more options
!               25 Feb 2003 (RRD) - ichar function replacment
!               11 Aug 2003 (RRD) - array limit test
!-------------------------------------------------------------------------------

      INTEGER, PARAMETER :: MAXPTS=300000

      REAL               :: RVAL(maxpts)
      CHARACTER(1)       :: buff(maxpts)
      CHARACTER(40)      :: fname

      WRITE(*,*)'Enter file name: '
      READ(*,'(a)')fname
      OPEN(10,file=fname,form='unformatted',access='direct',recl=1)
!###  OPEN(50,file='unpacker.txt')

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

!---->product definition "section 1" (unkown length - usually 28)
      krec=krec+8
      call fcread(10,buff,krec,28,kret)
      call sect1(buff,28,klen,kbms,kvar,kfact)

!---->grid description "section 2"
      if(kbms.eq.128.or.kbms.eq.192)then
         krec=krec+klen
         call fcread(10,buff,krec,32,kret)
         call sect2(buff,11,klen,nipt,njpt)
         nxyp=nipt*njpt
         if(nxyp.gt.maxpts)then
            write(*,*)'Data ',nxyp,'  exceeding compiled limit ',maxpts 
            stop
         end if
      end if

!---->bit map "section 3"
      if(kbms.eq.64.or.kbms.eq.192)then
         krec=krec+klen
         call fcread(10,buff,krec,6,kret)
         call sect3(buff,6,klen)
      end if

!---->binary data "section 4"
      krec=krec+klen
      call fcread(10,buff,krec,11,kret)
      call sect4(buff,11,klen,mfact,refval,nbpv)

      krec=krec+11
      call fcread(10,buff,krec,(klen-11),kret)

!###  if(kvar.eq.83)then

      CALL DECODE(buff,(klen-11),rval,nipt,njpt,nxyp,kfact,  &
                  mfact,refval,nbpv,cmin,cmax)

      write(*,'(I5,4(3X,A,E14.6))')                       &
         kvar,'Max:',cmax,'Min:',cmin,'V(1):',rval(1),'V(x):',rval(nxyp)
 
!###  end if

!---->grib end "section 5"
      krec=krec+(klen-11)
!     skip 4 octets
      krec=krec+4
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
      klen=jchar(buff(7))+jchar(buff(6))*256+jchar(buff(5))*65536
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT1(buff,nbyte,klen,kbms,kvar,kfact)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      kbms=jchar(buff(8))
      kvar=jchar(buff(9))
      kfact=jchar(buff(28))+jchar(buff(27))*256
      if(iand(kfact,32768).ne.0)kfact=-iand(kfact,32767)
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT2(buff,nbyte,klen,nipt,njpt)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      nipt=jchar(buff(8))+jchar(buff(7))*256
      njpt=jchar(buff(10))+jchar(buff(9))*256
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT3(buff,nbyte,klen)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE SECT4(buff,nbyte,klen,mfact,refval,nbpv)
      CHARACTER buff(nbyte)*1
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
      klen=jchar(buff(3))+jchar(buff(2))*256+jchar(buff(1))*65536

!     binary scale factor
      mfact=jchar(buff(6))+jchar(buff(5))*256
      if(iand(mfact,32768).ne.0)mfact=-iand(mfact,32767)

!     floating point value consists of 7 bit exponent and 24 bit mantissa
      mexp=iand(jchar(buff(7)),127)
      mant=jchar(buff(10))+jchar(buff(9))*256+jchar(buff(8))*65536

!     value equals 2^(-24)*MANT*16^(MEXP-64)
      refval=16.0**(mexp-64)*mant/16777216.0
      if(iand(jchar(buff(7)),128).ne.0)refval=-refval

!     number of bits per variable
      nbpv=jchar(buff(11))
      RETURN
      END

!-----------------------------------------------------------------

      SUBROUTINE DECODE(buff,nbyte,rval,nipt,njpt,nxyp,kfact,  &
                        mfact,refval,nbpv,cmin,cmax)

      CHARACTER buff(nbyte)*1
      REAL rval(nxyp)
      CHARACTER mychr*1
      INTEGER   jchar
      JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

      cmin= 1.0E+25
      cmax=-1.0E+25

      kv=0

      do jj=1,njpt
      do ii=1,nipt

         kv=kv+1
         rsum=0

         if(nbpv.gt.0)then
            ksum=0

!           compute starting (1) and ending (2) bit position
            kbit1=(kv-1)*nbpv+1
            kbit2=kv*nbpv

!           set the starting bit exponent as high value
            mexp=nbpv

!           loop through each bit left (high) to right (low)
            do kb=kbit1,kbit2

!              leftmost bit has value of 2^(nbpv-1)
               mexp=mexp-1

!              compute which byte contains the bit "kb"
               kbyte=(kb-1)/8+1

!              determine the relative bit position for "kb" in kbyte
               lbit=kb-(kbyte-1)*8

!              test if lbit in kbyte is turned on
               kbon=iand(jchar(buff(kbyte)),2**(8-lbit))

!              then that bit value to a final integer
               if(kbon.ne.0)ksum=ksum+2**mexp
            end do

!           compute scaled data point
            if(ksum.gt.0) rsum = float(ksum) * 2.0**mfact
         end if

!        real value
         rval(kv)= (refval+rsum) / 10.0**kfact

         cmax=max(cmax,rval(kv))
         cmin=min(cmin,rval(kv))

!###     write(50,'(2I4,F6.2)')II,JJ,RVAL(KV) 

      end do
      end do

      RETURN
      END
