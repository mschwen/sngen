***************************************************************************
* Filename: en_nu_d_n_cc.f
***************************************************************************
*
*  This routine produces both randomly sampled neutrons for the reaction:
*  nu_e_bar + d --> n + n + e+
*
*  This routine relies on the fact that the positron has already been
*  generated (in this case, it is done in 'en_nu_d_cc.f')
*  It also assumes that the di-neutron system left over after the CC
*  interaction decays isotropically in its frame
*
*  Inputs:
*	o  neutrino energy (MeV)
*  Outputs:
*	o  sampled final state neutron kinetic energy (MeV)
*          (either neutron #1 or neutron #2)
*
***************************************************************************

      SUBROUTINE en_nu_d_n_cc(E_nu,E_part,redo)
   
      IMPLICIT none

      INCLUDE 'physics_constants.inc'

C   Common Block Declarations

      COMMON /nu_d_positron_stuff/ E_e, sn_dir_cos, part_dir 
      REAL*8 E_e, sn_dir_cos(3), part_dir(3)

c      COMMON /nu_d_dir_e/ cos_theta_e
c      REAL*8 cos_theta_e

      REAL*8 cos_theta_n1                   ! needed by 'dir_nu_d_n_nc.f'
      COMMON /nu_d_dir_n1/ cos_theta_n1

      DOUBLE PRECISION En1, En2             ! needed by 'dir_nu_d_n_nc.f'
      COMMON /nu_d_E_neutrons/ En1, En2

C   Global Variable Declarations

      REAL*8  E_nu, E_part
      INTEGER redo

C   Local Variable Declarations

      DOUBLE PRECISION dir_nu(3), dir_e(3)

      DOUBLE PRECISION p_e
      DOUBLE PRECISION E_dn, p_dn(3), p_dn2, m_dn, m_dn2
c      DOUBLE PRECISION rtemp
c      DOUBLE PRECISION E_dn, p_dn(3), p_dn2, p_dn2_new, m_dn, m_dn2
      DOUBLE PRECISION p_n1_cms, e_n1_cms

      DOUBLE PRECISION ran_val, cos_ran, cos_cms

      DOUBLE PRECISION pz_n1_cms

      DOUBLE PRECISION gamma, beta, pz_n1_lab

      INTEGER en_neutron_count
      DOUBLE PRECISION E_n1, E_n2

      DOUBLE PRECISION E_diff

      SAVE en_neutron_count
      SAVE E_n1,E_n2               ! save for neutron #2 dir kinematics
      DATA en_neutron_count /0/
      DATA E_n1, E_n2 /2 * 0.d0/

***************************************************************************

      redo = 0

****************************************************************************
*  Neutron #1 (kinetic) energy
*  [for nu_e_bar + d --> n + n + e+  CC reaction]
****************************************************************************

      IF (en_neutron_count .EQ. 0) THEN      

*** supernova direction cosines (in SNO lab frame)

         dir_nu(1) = sn_dir_cos(1)
         dir_nu(2) = sn_dir_cos(2)
         dir_nu(3) = sn_dir_cos(3)

*** positron direction cosines (in lab frame)

         dir_e(1) = part_dir(1)
         dir_e(2) = part_dir(2)
         dir_e(3) = part_dir(3)

****************************************************************************
*  First calculate di-neutron invariant mass, using conservation of
*  energy and momentum in the lab frame
****************************************************************************

***  positron momentum (lab frame)
 
         p_e = dsqrt(E_e**2 - m_electron**2)

***  di-neutron energy (using conservation of total energy in lab frame)

         E_dn = E_nu + m_deuteron - E_e

         E_diff = E_nu - (E_e - m_electron)

***  di-neutron momenta (lab frame)

         p_dn(1) = dir_nu(1) * E_nu - dir_e(1) * p_e
         p_dn(2) = dir_nu(2) * E_nu - dir_e(2) * p_e
         p_dn(3) = dir_nu(3) * E_nu - dir_e(3) * p_e

         p_dn2 = p_dn(1)**2 + p_dn(2)**2 + p_dn(3)**2

c         the following gives the same result:
c         p_dn2_new = E_nu**2 + p_e**2 - 2.d0 * E_nu * p_e * cos_theta_e

***  di-neutron mass

         m_dn = dsqrt(E_dn**2 - p_dn2)
   
*** small fudge factor to save things from dying on 'p_n1_cms' 
*** for now until I get it fixed

           IF (m_dn .LT. (2.d0 * m_neutron)) THEN
              redo = 1
              RETURN
           ENDIF
    
         m_dn2 = m_dn * m_dn

**************************************************************************** 
*  Now let di-neutron decay into two neutrons in its frame (CMS)
****************************************************************************

*** Handle rounding errors

c         rtemp = max(m_dn2 - (2.d0*m_neutron)**2,m_electron**2)

*** Proceed with calculation

c         p_n1_cms = dsqrt(rtemp * m_dn2 / (2.d0 * dsqrt(m_dn2)))
         p_n1_cms = 0.5d0 * dsqrt(m_dn2 - (2.d0 * m_neutron)**2)
         E_n1_cms = dsqrt(p_n1_cms**2 + m_neutron**2)

****************************************************************************
*  Pick decay angle isotropically in di-neutron frame
*  z axis is along direction of di-neutron momentum
****************************************************************************

         cos_ran = ran_val(-1.d0,1.d0)
         cos_cms = cos_ran 

         pz_n1_cms = p_n1_cms * cos_cms

*****************************************************************************
*  Boost along di-neutron momentum---E_n1 is neutron energy in lab,
*  'cos_theta' is angle in lab relative to di-neutron direction.
*****************************************************************************

         gamma = E_dn / m_dn
         beta = dsqrt(gamma**2 - 1.d0) / gamma
         E_n1 = gamma * (beta * pz_n1_cms + E_n1_cms)
         pz_n1_lab = gamma * (beta * E_n1_cms + pz_n1_cms)

         En1 = E_n1   ! fill COMMON block

* Calculate neutron #1 (kinetic) energy (in MeV):

         E_part = E_n1 - m_neutron

* Calculate neutron #1 cos_theta scattering direction
* (used for direction calculation in 'dir_nu_d_n_cc.f')
 
         cos_theta_n1 = pz_n1_lab / dsqrt(E_n1**2 - m_neutron**2) 
c         cos_theta_n1 = pz_n1_lab / (E_n1**2 - m_neutron**2) ! error!

         IF (abs(cos_theta_n1) .GT. 1.d0) cos_theta_n1 = -1.d0 ! ensure physical 
* if 'cos_theta' is outside the bounds, set it to +/- 1 
* (depending on what* the sign of the original value was)
* Note: this implementation was responsible for excess events at cos(theta)~0
c         IF (abs(cos_theta_n1) .GT. 1.d0) THEN
c            cos_theta_n1 = sign(1.d0,cos_theta_n1)
c         ENDIF

* Increment the neutron counter for the next neutron in this neutrino event

         en_neutron_count = 1

****************************************************************************
*  Neutron #2 (kinetic) energy
*  [for nu_e_bar + d --> n + n + e+  CC reaction]
****************************************************************************

      ELSEIF (en_neutron_count .EQ. 1) THEN

         E_n2 = E_nu + m_deuteron - E_e - E_n1

***  Fill common block

         En2 = E_n2

***  Calculate neutron #2 (kinetic) energy (in MeV):

         E_part = E_n2 - m_neutron

***  Handle rounding errors (see notes in SNOMAN en_cc_neutron.for)

c         E_part = max(E_part,m_electron)

***  Reset 'en_neutron_counter' back to '0' for the next neutrino event

         en_neutron_count = 0

      ENDIF

      RETURN
      END
