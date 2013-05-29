***************************************************************************
* Filename: dir_nu_p_es.f
***************************************************************************
*
*   Inputs:
*	o  total neutrino energy (MeV)
*	o  total electron energy  (MeV)
*	o  supernova direction cosines
*   Outputs:
*	o  sampled final state particle directional cosines
*
***************************************************************************

      SUBROUTINE dir_nu_p_es(E_nu,E_proton,sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'

C   Global Variable Declarations

      REAL*8 E_nu
      DOUBLE PRECISION E_proton
      REAL*8 sn_dir_cos(*)
      REAL*8 part_dir(3)

C   Local Variable Declarations

      DOUBLE PRECISION T_p

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION cos_theta_p, theta_p, azimuth_ran, azimuth

***************************************************************************

***  Calculate the proton  kinetic energy

      T_p = E_proton

***  Calculate the scattering angle, 'theta', from the kinematics

      cos_theta_p = dsqrt(T_p / (T_p + 2.d0 * m_proton)) 
     +              * ((E_nu + m_proton) / E_nu)

      IF (cos_theta_p .GT. 1.d0) THEN
         WRITE (*,*) '*** dir_nu_p_es.f: Something went wrong!'
         WRITE (*,*) '*** dir_nu_p_es.f: cos_theta_p =', cos_theta_p
      ENDIF

      theta_p = dacos(cos_theta_p)

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi.

      azimuth_ran = ran_val(0.d0,1.d0)
      azimuth = (2.d0 * pi) * azimuth_ran

***************************************************************************
*  Combine supernova angle with electron scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta_p,azimuth,part_dir)

***************************************************************************

      RETURN
      END
