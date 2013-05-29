***************************************************************************
* Filename: dir_nu_d_nc.f
***************************************************************************
*
*   Inputs:
*	o  supernova direction cosines
*   Outputs:
*	o  sampled final state particle directional cosines
*
***************************************************************************

      SUBROUTINE dir_nu_d_nc(sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C   Common Block Declarations

      REAL*8 cos_theta                  ! comes from 'en_nu_d_nc.f'
      COMMON /nu_d_dir_stuff/ cos_theta

C   Global Variable Declarations

      REAL*8 sn_dir_cos(*)
      REAL*8 part_dir(3)

C   Local Variable Declarations

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION theta, azimuth_ran, azimuth

*****************************************************************************

      theta = dacos(cos_theta)

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi.

      azimuth_ran = ran_val(0.d0,1.d0)
      azimuth = (2.d0 * pi) * azimuth_ran

***************************************************************************
*  Combine supernova angle with the final state particle scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta,azimuth,part_dir)

***************************************************************************

      RETURN
      END
