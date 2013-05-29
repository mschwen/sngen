**********************************************************************
*  Filename: equator2horizon.f
***********************************************************************
*
*  Transforms the "equatorial" coordinate system to the "local horizon"
*  coordinate system of SNO. 
*
***********************************************************************

      SUBROUTINE equator2horizon(hourang,declin,det_latitude,
     +                           altitude,azimuth)

      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C  Global Variable Declarations (in)

      DOUBLE PRECISION hourang, declin, det_latitude, altitude

C  Global Variable Declarations (out)

      DOUBLE PRECISION azimuth

C  Local Variable Declarations

      DOUBLE PRECISION over_piby2, under_piby2      
      DOUBLE PRECISION abs_altitude,
     +                 sin_azimuth, abs_sin_azimuth
      DOUBLE PRECISION lhs, rhs, prod

*************************************************************************

C       ANGLES ALL GIVEN IN RADIANS
C       ALTITU  :  Horizon altitude
C       AZIMUT  :  Horizon Azimuth
C       LATITU  :  DETECTOR LATITUDE
C       HOURAN  :  Hour Angle in radians
C       DECLIN  :  Declination

*************************************************************************

c      WRITE (*,*) 'Converting to local horizon coordinates...'

**************************************************************************
*   Find the horizon angles (altitude and azimuth) in the SNO detector
*   frame. 
*
*   This routine will run into trouble if "cos(altitude)=piby2". What that
*   means is that the supernova is directly on the zenith or nadir of 
*   the detector. Need to specify that information in the code. 
**************************************************************************

* Bounds for REAL*8 number around 90 degrees

      over_piby2  = 90.0000001 * deg2rad
      under_piby2 = 89.9999999 * deg2rad

* Conversion formulae

      altitude =  asin(sin(det_latitude) * sin(declin) 
     +          + cos(det_latitude) * cos(declin) * cos(hourang))

      abs_altitude = abs(altitude)

      IF (      abs_altitude .GT. under_piby2
     +    .AND. abs_altitude .LT. over_piby2 ) THEN

         azimuth = piby2

      ELSE
         sin_azimuth = -cos(declin) * sin(hourang) / cos(altitude)

         abs_sin_azimuth = abs(sin_azimuth)

         IF (abs_sin_azimuth .LT. 1.d0) THEN
            azimuth = asin(sin_azimuth)
         ELSE
            azimuth = sign(piby2,sin_azimuth)  ! bug fix Dec 8, 2001
*            azimuth = piby2
         ENDIF
      ENDIF
 
*  Azimuth: -PI/2 ----> + PI/2          HA : -PI/2 ----> + PI/2
*  REMOVE THE AMBIGUITY ARISING FROM ARCSIN
 
      lhs = cos(altitude) * cos(azimuth)
      rhs = sin(declin) * cos(det_latitude)
     +     - cos(declin) * cos(hourang) * sin(det_latitude)
 
      IF (lhs .EQ. 0.d0) THEN
        lhs = lhs + 1.d0
      ENDIF
 
      prod = lhs * rhs

      IF (prod .NE. (abs(prod))) THEN
        IF (azimuth .GT. 0.d0) THEN
          azimuth = pi - azimuth
        ELSE
          azimuth = -pi - azimuth
        ENDIF
      ENDIF
 
      IF (azimuth .LT. 0.d0) THEN
        azimuth  = azimuth + 2.d0 * pi
      ENDIF

c      WRITE (*,*) 'altitude, azimuth =', altitude, azimuth

      RETURN
      END
