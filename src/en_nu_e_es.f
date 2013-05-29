***************************************************************************
* Filename: en_nu_e_es.f
***************************************************************************
*
*   Inputs:
*	o  interaction id number (not array)
*	o  delayed neutrino interaction time (ms)
*	o  neutrino energy (MeV)
*   Outputs:
*	o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_e_es(interaction,E_nu,E_part)
   
      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'physics_constants.inc'

C   Global Variable Declarations

      INTEGER interaction
      REAL*8 E_nu
      DOUBLE PRECISION E_part

C   Local Variable Declarations

      INTEGER nu_type

      REAL*8 g_L, g_R, g_L2, g_R2

      DOUBLE PRECISION get_esdif_rad
      DOUBLE PRECISION T_e_max, T_e

      LOGICAL done

      DOUBLE PRECISION ran_val, xsect_ran

      REAL*8 xsect_max, xsect_diff, xsect_ratio

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
*  Similar reactions in LAB and 208Pb
************************************************************************

      CALL get_nu_type(interaction,nu_type)

***************************************************************************
*  Calculate d_(xsection)/d_T for electron recoil energy, e_electron
***************************************************************************

* radiative corrections
      IF (     nu_type .EQ. id_nu_e
     +    .OR. nu_type .EQ. id_nu_mu
     +    .OR. nu_type .EQ. id_nu_tau) THEN

* nu_e
c         g_L = (+0.5d0 + sin2theta)
c         g_R = sin2theta

* nu_mu, nu_tau
c         g_L = (-0.5d0 + sin2theta)
c         g_R = sin2theta

         T_e  = 0.00001d0
         xsect_max = get_esdif_rad(nu_type, E_nu, T_e)

      ELSE
         IF (nu_type .EQ. id_nu_ebar) THEN

            g_L = sin2theta
            g_R = (+0.5d0 + sin2theta)

         ELSEIF (     nu_type .EQ. id_nu_mubar
     +           .OR. nu_type .EQ. id_nu_taubar) THEN

            g_L = sin2theta
            g_R = (-0.5d0 + sin2theta)

         ENDIF

         g_L2 = g_L * g_L
         g_R2 = g_R * g_R

**** maximum differential cross section (without constant sig_e factor)

         xsect_max = g_L2 + g_R2

      ENDIF  ! radiative corrections

***************************************************************************
*  Calculate differential cross section based on random cross section and
*  decide whether to keep cross section evaluated at recoil energy
***************************************************************************

***  Calculate maximum allowed electron kinetic energy

      T_e_max = 2.d0 * E_nu*E_nu / (2.d0 * E_nu + m_electron)

      done = .false.
         
      DOWHILE (.NOT. done)

         xsect_ran = ran_val(0.d0,1.d0)
         T_e = T_e_max * xsect_ran

* calculate test cross section

         IF (     nu_type .EQ. id_nu_e
     +          .OR. nu_type .EQ. id_nu_mu
     +          .OR. nu_type .EQ. id_nu_tau) THEN

* radiative corrections
            xsect_diff = get_esdif_rad(nu_type, E_nu, T_e)

         ELSE

* no radiative corrections
            xsect_diff =  g_L2 + g_R2 * (1.d0 - T_e / E_nu)**2 
     +                  - g_L * g_R * m_electron * T_e / E_nu**2

         ENDIF  ! radiative corrections

         xsect_ratio = xsect_diff / xsect_max
     
         xsect_ran = ran_val(0.d0,1.d0)
         IF (xsect_ran .GT. xsect_ratio) THEN
            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO

*** Calculate electron kinetic energy

      E_part = T_e + m_electron

C      xsect_diff = sig_e * xsect_diff
C      WRITE (*,*) E_part

***************************************************************************

      RETURN
      END
