***********************************************************************
*  Filename: sn_direction.f
***********************************************************************
*
*  Converts the "equatorial" coordinate system to the "local horizon"  
*  coordinate system of SNO. After that the direction cosines in       
*  SNO's frame is calculated. Equations are taken from Spherical       
*  Astronomy by Woolard and Clemence or Lang
*
***********************************************************************

      SUBROUTINE sn_direction(sn_dir_cos)

      IMPLICIT none

      INCLUDE 'sno_coordinates.inc'
      INCLUDE 'trig_constants.inc'

C  Common Block Declarations

      COMMON /direction/ deg_declin, deg_rascen,
     +                   dsec, imin, ihr, iday, imonth, iyear
      DOUBLE PRECISION   deg_declin, deg_rascen
      DOUBLE PRECISION   dsec
      INTEGER            imin, ihr, iday, imonth, iyear

C  Global Variable Declarations (out)

      DOUBLE PRECISION   sn_dir_cos(3)

C  Local Variable Declarations

      DOUBLE PRECISION declin, rascen, hourang,
     +                 det_latitude, det_longitude,
     +                 altitude, azimuth
      DOUBLE PRECISION dusec
      DOUBLE PRECISION nday, gmst
      INTEGER          idst, ideltas

**************************************************************************

*   Convert "celestial" coordinate system angles from degrees to radians.

      declin = deg_declin * deg2rad
      rascen = deg_rascen * deg2rad

*   Convert SNO detector coordinates from degrees to radians.

      det_latitude  = SNO_latitude * deg2rad
      det_longitude = SNO_longitude * deg2rad

***************************************************************************
*   Convert Sudbury time to universal time
***************************************************************************

* decide if we need to correct for daylight savings time
* (-1 = EDT  (GMT minus 4 hours); 0 = EST (GMT minus 5 hours))

      CALL msw_ee_calc_dst(iyear,imonth,iday,ihr,idst)

      ideltas = 3600 * (5 + idst)

* get Greenwich mean standard time

      dusec =  dble(ihr) * 3600.d0             ! number of seconds in day
     +       + dble(imin) * 60.d0 
     +       + dsec
     +       + dble(ideltas)

      IF (int(dusec) .GT. 86400) THEN
         iday = iday + 1
         dusec = dusec - 86400.d0
      ENDIF

c      WRITE (*,*) 'iyear, imonth, iday, dusec =',
c     +             iyear, imonth, iday, dusec

      CALL daynum(iyear,imonth,iday,nday)      ! get day number (SNOMAN)
      CALL utgmt(iyear,imonth,nday,dusec,gmst) ! get gmst (SNOMAN)

****************************************************************************
*   Convert from "celestial" coordinates to "equatorial" coordinates
****************************************************************************

      CALL celest2equator(rascen,gmst,det_longitude,hourang)

****************************************************************************
*   Convert from "equatorial" coordinates to "local horizon" coordinates
****************************************************************************

      CALL equator2horizon(hourang,declin,det_latitude,altitude,azimuth)

****************************************************************************
*   Determine the SNO detector coordinate system directional cosines
****************************************************************************

      CALL horizon2det(altitude,azimuth,sn_dir_cos)

      RETURN
      END
