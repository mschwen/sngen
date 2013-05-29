***************************************************************************
* Filename: dir_nu_e_es.f
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

      SUBROUTINE dir_nu_e_es(E_nu,e_electron,sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'

C   Global Variable Declarations

      REAL*8 E_nu
      DOUBLE PRECISION E_electron
      REAL*8 sn_dir_cos(*)
      REAL*8 part_dir(3)

C   Local Variable Declarations

      DOUBLE PRECISION T_e

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION cos_theta_e, theta_e, azimuth_ran, azimuth

***************************************************************************

************************************************************************
*   1 nu_e + e- -> nu_e + e-  H2O (ES)
*   2 nu_e + e- -> nu_e + e-  D2O (ES)
*  10 nu_e_bar + e- -> nu_e_bar + e-  H2O (ES)
*  11 nu_e_bar + e- -> nu_e_bar + e-  D2O (ES)
*  21 nu_mu + e- -> nu_mu + e-  H2O (ES)
*  22 nu_mu + e- -> nu_mu + e-  D2O (ES)
*  30 nu_mu_bar + e- -> nu_mu_bar + e-  H2O (ES)
*  31 nu_mu_bar + e- -> nu_mu_bar + e-  D2O (ES)
*  39 nu_tau + e- -> nu_tau + e-  H2O (ES)
*  40 nu_tau + e- -> nu_tau + e-  D2O (ES)
*  48 nu_tau_bar + e- -> nu_tau_bar + e-  H2O (ES)
*  49 nu_tau_bar + e- -> nu_tau_bar + e-  D2O (ES)
************************************************************************

***  Calculate the electron kinetic energy

      T_e = E_electron - m_electron

***  Calculate the scattering angle, 'theta', from the kinematics
***  (see for example Bahcall, Beacom & Vogel astro-ph/9811350v2)

      cos_theta_e = dsqrt(T_e / (T_e + 2.d0 * m_electron)) 
     +              * ((E_nu + m_electron) / E_nu)

      IF (cos_theta_e .GT. 1.d0) THEN
         WRITE (*,*) '*** dir_nu_e_es.f: Something went wrong!'
         WRITE (*,*) '*** dir_nu_e_es.f: cos_theta_e =', cos_theta_e
      ENDIF

      theta_e = dacos(cos_theta_e)

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi.

      azimuth_ran = ran_val(0.d0,1.d0)
      azimuth = (2.d0 * pi) * azimuth_ran

***************************************************************************
*  Combine supernova angle with electron scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta_e,azimuth,part_dir)

***************************************************************************

      RETURN
      END
