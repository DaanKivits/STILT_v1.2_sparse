      include 'cmapf.fi'
      REAL STCPRM(k_maparam)
      CALL STLMBR (STCPRM, EQVLAT(60.,40.), - 70.)
C	CALL STCM2P (STCPRM, 1.,1., 45.,-87., 10.,1., 45., - 85.)
      CALL STCM1P (STCPRM, 1.,1., 45.,-87., 45.,-70., 17.5, 0.)
c***********
      call stlmbr (STCPRM,30., -80.)
      call stcm2p (STCPRM, 1.,1., 24.,-86., 11.,6., 42.,-71.5)
      write(*,*) 'Gamma =',STCPRM(k_gama)
      write(*,*) 'grid size at 30 is', cgszll (STCPRM, 30.,-80.)
c********
      do lat = 0,90,10
         do long=-87,-85
            write(*,*) ' -------------------------------------'
            xlat = lat
            xlong = long
            call cll2xy(STCPRM, xlat,xlong, x,y)
            write (*,*) 'from lat,long=',xlat,',',xlong
            write (*,*) '       to x,y=',x,y
            call cxy2ll(STCPRM, x,y, xlat1,xlong1)
            write (*,*) '   from x,y=',x,',',y
            write (*,*) 'to lat,long=',xlat1,xlong1
            write(*,*) '           ---'
            ue = 30.
            vn = 40.
            call cc2gxy(STCPRM, x,y, ue,vn, ug,vg)
            write(*,*) 'xy- from compass wind',ue,',',vn
            write(*,*) '    to grid wind', ug,',',vg,' (',
     c	    180. - dgprad*atan2(ug,-vg),' deg)'
            call cc2gll(STCPRM, xlat,xlong, ue,vn, ug,vg)
            write(*,*) 'll- from compass wind',ue,',',vn
            write(*,*) '    to grid wind', ug,',',vg,' (',
     c	    180. - dgprad*atan2(ug,-vg),' deg)'
            call cg2cxy(STCPRM, x,y, ug,vg, ue,vn)
            write(*,*) 'xy- from grid wind',ug,',',vg
            write(*,*) '    to compass wind',ue,',',vn,' (',
     c      180. - dgprad*atan2(ue,-vn),' deg)'
            call cg2cll(STCPRM, xlat,xlong, ug,vg, ue,vn)
            write(*,*) 'll- from grid wind',ug,',',vg
            write(*,*) '    to compass wind',ue,',',vn,' (',
     c      180. - dgprad*atan2(ue,-vn),' deg)'
            write(*,*) 'xy- gridsize =',cgszxy(STCPRM, x,y),' km'
            write(*,*) 'll- gridsize =',cgszll(STCPRM, xlat,xlong),
     c      ' km'
            call cpolxy(STCPRM, x,y, enx,eny,enz)
            write (*,*) 'xy- polar components=',enx,',',eny,',',enz
            call cpolll (STCPRM, xlat,xlong, enx,eny,enz)
            write (*,*) 'll- polar components=',enx,',',eny,',',enz
            call ccrvxy(STCPRM, x,y, gx,gy)
            write(*,*) 'xy- curvature vector',gx,',',gy
            call ccrvll(STCPRM, xlat,xlong, gx,gy)
            write(*,*) 'll- curvature vector',gx,',',gy
         enddo
      enddo
      STOP
      END
