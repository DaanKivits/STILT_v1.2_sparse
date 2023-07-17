      function eqvlat(lat1,lat2)
C*  Written 12/21/94 by Dr. Albion Taylor
C*
C*    This function is provided to assist in finding the tangent latitude
C*    equivalent to the 2-reference latitude specification in the legend
C*    of most lambert conformal maps.  If the map specifies "scale
C*    1:xxxxx true at 40N and 60N", then eqvlat(40.,60.) will return the
C*    equivalent tangent latitude.
C*  INPUTS:
C*    lat1, lat2:  The two latitudes specified in the map legend
C*  RETURNS:
C*    the equivalent tangent latitude
C*  EXAMPLE:  stcmap(& strcmp, eqvlat(40.,60.), 90.)
C*/
C/*   Changes Made May 9, 2003 to accomodate the following special
C*   situations:
C*   1. if lat1 == lat2, returned value will be lat1 (reduced to between -90.
C*      and 90.).
C*   2. If either lat1 or lat2 is 90. (or -90.) then 90. (or -90.) will be
C*      returned.  This reflects the fact that, for y fixed (-90. < y < 90.),
C*      as x -> 90. (or -90.), eqvlat(x,y) ->90. (or -90.)  This limiting
C*      case of tangent latitude 90. is a polar stereographic projection,
C*      for which the scale at 90. is a maximum, and therefore greater than the
C*      other latitude y.  Thus, eqvlat(90.,60.) returns 90., although the
C*      scale at 90. will be greater than at 60. For eqvlat(90.,-90.), the
C*      limit process is ambiguous; for the sake of symmetry, such a case
C*      will return 0.
C*/
      include 'cmapf.fi'

	real lat1,lat2
      parameter (fsm = 1.e-3)
        slat1 = sin(radpdg * lat1)
        slat2 = sin(radpdg * lat2)
C reorder slat1, slat2
        if (  slat1 .lt. slat2) then
          temp = slat1
          slat1 = slat2
          slat2 = temp
        endif
C/*  Take care of special cases first */
        if (slat1 .eq. slat2) then
          eqvlat = asin(slat1) * DGPRAD
          return
        endif
        if (slat1 .eq. -slat2 ) then
          eqvlat = 0.
          return
        endif
        if (slat1 .ge. 1.) then
          eqvlat = 90.
          return
        endif
        if (slat2 .le. -1.) then
          eqvlat = -90.
          return
        endif
C/* Compute al1 = log((1. - slat1)/(1. - slat2))/(slat1 - slat2) */
        tau = (slat1 - slat2)/(2. - slat1 - slat2)
        if ( tau .gt. FSM )then
          al1 = log((1. - slat1)/(1. - slat2)) / (slat1 - slat2)
        else
          tau = tau * tau
          al1 = -2./(2. - slat1 - slat2) *
     a                  (1.    + tau *
     b                  (1./3. + tau *
     c                  (1./5. + tau *
     d                  (1./7. + tau))))
        endif
C/* Compute al2 = log((1. + slat1)/(1. + slat2))/(slat1 - slat2) */
        tau = (slat1 - slat2)/(2. + slat1 + slat2)
        if ( tau .gt. FSM ) then
          al2 = log((1. + slat1)/(1. + slat2)) / (slat1 - slat2)
        else
          tau = tau * tau
          al2 =  2./(2. + slat1 + slat2) *
     a                  (1.    + tau *
     b                  (1./3. + tau *
     c                  (1./5. + tau *
     d                  (1./7. + tau))))
        endif
        eqvlat = asin ((al1 + al2) / (al1 - al2)) * DGPRAD
        return
      end

