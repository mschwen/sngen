**********************************************************************
*  Filename: horizon2det.f
**********************************************************************
*
* Converts the zenith and azimuth in the horizon system to supernova  
* direction cosines (directional vector) with respect to the centre   
* of the SNO detector (in PSUP coordinates).
* 
* Use SNOMAN direction references from "solar.for" and correction from
* SNO-STR-99-004
*
**********************************************************************

      SUBROUTINE horizon2det(altitude,azimuth,dir_cos)

      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C  Global Variable Declarations (in)

      DOUBLE PRECISION altitude, azimuth

C  Global Variable Declarations (out)

      DOUBLE PRECISION dir_cos(3)

C  Local Variable Declarations

      DOUBLE PRECISION theta, phi
      INTEGER i

****************************************************************************
*
*     dir_cos(1) = directional cosine along x-axis
*     dir_cos(2) = directional cosine along y-axis
*     dir_cos(3) = directional cosine along z-axis
*
*     pi = your favourite piece of number
*     piby2 = pi divided by 2
*
***************************************************************************

* Calculate the supernova direction cosines along each of the 3 directions:
*     PSUP North is 49.58 deg West of True (Geographic) North
*     PSUP East is  40.42 deg East of True North
*     'x' --> points East
*     'y' --> points to true North (GEOGRAPHIC NORTH)
*     'z' --> points up toward zenith (PARALLEL TO DETECTOR AXIS)
*     x = y = z = 0 --> centre of the acrylic vessel
*     'theta' --> 'zenith' angle measured from the z-axis (0 to pi)
*     'phi'   --> 'azimuthal' angle measured eastwards from the y-axis 
*                  (0 to 2*pi)

* Convert altitude angle to zentih angle

      theta = piby2 - altitude

* Convert azimuthal angle to PSUP coordinates

      phi = azimuth + 0.86533d0

c      WRITE (*,*) 'theta, phi =', theta, phi

* Calculate direction cosines

      dir_cos(1) = sin(theta) * sin(phi)
      dir_cos(2) = sin(theta) * cos(phi)
      dir_cos(3) = cos(theta)

* Calculate the neutrino direction cosines along each of the 3 directions
* (the negative accounts for the fact that the neutrino momentum direction
*  is opposite to the supernova)

      DO 10 i = 1,3
         dir_cos(i) = -dir_cos(i)
c         write (*,*) 'Toward SN: dircos(',i,') =', -dir_cos(i)! output SN dircos
10    CONTINUE

      RETURN
      END
