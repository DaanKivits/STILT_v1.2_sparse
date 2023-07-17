      DIMENSION stdgrd(11)
      CALL STLMBR(STDGRD, 90., -80.)
      CALL STCM1P(STDGRD, 33.,33., 90.,-80., 60.,-80., 381.,0.)
      write(*,* ) 'Initialization complete'
      DELTAX=1.
      DO K=1,26
        CALL CXY2LL(STDGRD, 33.-DELTAX,33.-DELTAX,XLAT,XLONG)
        WRITE(*,*) DELTAX, XLAT, XLONG, 
     C       DELTAX*CGSZXY(STdgrd, 33.-DELTAX,33.-DELTAX)
        IF (K.EQ.25) DELTAX=0.
        DELTAX=.5*DELTAX
      ENDDO
      END
