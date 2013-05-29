**********************************************************************
*  Filename: celest2equator.f
***********************************************************************
*
*  Given a time (GMST), this routine transforms the "celestial" coordinate
*  system (ra,dec) to the "equatorial" coordinate system (hour angle,dec).  
*
*  Note: All angles are in radians.
*
***********************************************************************

      SUBROUTINE celest2equator(rascen,gmst,det_longitude,hourang)

      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C  Global Variable Declarations

      DOUBLE PRECISION  rascen, det_longitude, hourang
      DOUBLE PRECISION  gmst

C  Local Variable Declarations
      
      DOUBLE PRECISION  lst
      DOUBLE PRECISION  dayby2         ! number of seconds in 12 hours
      PARAMETER         (dayby2 = 43200.0d+0)

*************************************************************************

C       ANGLES ALL GIVEN IN RADIANS
C       HOURAN  :  Hour Angle in radians
C       RASCEN  :  Right Ascension.
C       LONGIT  :  Detector longitude.
C       GMST    :  Greenwich Mean Sidereal Time.
C       LST     :  Local Mean Sidereal Time in seconds.

**************************************************************************

c      WRITE (*,*) 'Converting to equatorial coordinates...'

* Calculate local mean sidereal time (in seconds)

      lst       = gmst - (det_longitude / pi) * dayby2

* Calculate hour angle (in radians)

      hourang   = (lst / dayby2) * pi - rascen

* Assign hour angle values to the other quadrants beyond pi

      DOWHILE (hourang .LE. (-pi))
         hourang = hourang + 2.d0*pi
      END DO

      DOWHILE (hourang .GT. (pi))
         hourang = hourang - 2.d0*pi
      END DO

c      WRITE (*,*) 'gmst, lst, hourang =', gmst, lst, hourang

      RETURN
      END
