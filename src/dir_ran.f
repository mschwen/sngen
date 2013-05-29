***************************************************************************
* Filename: dir_ran.f
***************************************************************************
*
*   Inputs:
*	o  supernova direction cosines
*   Outputs:
*	o  sampled final state particle directional cosines
*
***************************************************************************

      SUBROUTINE dir_ran(sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C   Global Variable Declarations

      REAL*8 sn_dir_cos(*)
      REAL*8 part_dir(3)

C   Local Variable Declarations

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION theta, azimuth_ran, azimuth, theta_ran
      INTEGER i      

*****************************************************************************

***  Choose the axial scattering angle uniformly between 0 and 2*pi.

      theta_ran = ran_val(0.d0, 1.d0)
      theta = (2.d0 * pi) * theta_ran

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi.

      azimuth_ran = ran_val(0.d0,1.d0)
      azimuth = (2.d0 * pi) * azimuth_ran

* This was used to create a file to test the randomness of ran_val; it should commented
C      OPEN (unit = 1, file =  'Random.out', status = 'new')
C        DO i = 1, 1000000
C           theta_ran = ran_val(0.d0, 1.d0)
C           theta = (2.d0 * pi) * theta_ran
C           WRITE (1,45) theta
C45         FORMAT (F14.1)
C        END DO
C      CLOSE (unit = 1)

***************************************************************************
*  Combine supernova angle with the final state particle scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta,azimuth,part_dir)

***************************************************************************

      RETURN
      END
