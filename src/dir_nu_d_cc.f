***************************************************************************
* Filename: dir_nu_d_cc.f
***************************************************************************
*
*   Inputs:
*	o  interaction id number (not array)
*	o  total neutrino energy (MeV)
*	o  total electron energy  (MeV)
*	o  supernova direction cosines
*   Outputs:
*	o  sampled final state particle directional cosines
*
***************************************************************************

      SUBROUTINE dir_nu_d_cc(E_part,sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'trig_constants.inc'

C   Common Block Declarations

      REAL*8 cos_theta_e
      COMMON /nu_d_dir_e/ cos_theta_e

***   Fill common blocks with values needed for neutron kinematics
      DOUBLE PRECISION nu_dir_cos(3), part_e_dir(3)
      COMMON /nu_d_positron_stuff/ E_e, nu_dir_cos, part_e_dir 

C   Global Variable Declarations

      REAL*8 E_part
      REAL*8 sn_dir_cos(*)
      REAL*8 part_dir(3)

C   Local Variable Declarations

      DOUBLE PRECISION E_e

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION theta_e, azimuth_ran, azimuth

      INTEGER ith_dir

****************************************************************************
*  Fill common block parameters needed for neutron calculations
****************************************************************************

      E_e = E_part

      DO 10 ith_dir = 1, 3
         nu_dir_cos(ith_dir) = sn_dir_cos(ith_dir)
         part_e_dir(ith_dir) = part_dir(ith_dir)
10    CONTINUE

*****************************************************************************

      theta_e = dacos(cos_theta_e)

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi

      azimuth_ran = ran_val(0.d0,1.d0)
      azimuth = (2.d0 * pi) * azimuth_ran

***************************************************************************
*  Combine supernova angle with the final state particle scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta_e,azimuth,part_dir)

***  fill array for common block for nu + d --> n + n + e+ reaction
      DO 20 ith_dir = 1, 3
         part_e_dir(ith_dir) = part_dir(ith_dir)
20    CONTINUE

***************************************************************************

      RETURN
      END
