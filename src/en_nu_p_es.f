***************************************************************************
* Filename: en_nu_p_es.f
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

      SUBROUTINE en_nu_p_es(interaction,E_nu,E_part)
   
      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'physics_constants.inc'

C   Global Variable Declarations

      INTEGER interaction
      REAL*8 E_nu
      DOUBLE PRECISION E_part

C   Local Variable Declarations

      INTEGER nu_type

      REAL*8 ca, cv

      DOUBLE PRECISION T_p_max, T_p

      LOGICAL done

      DOUBLE PRECISION ran_val, xsect_ran

      REAL*8 xsect_max, xsect_diff, xsect_ratio

***************************************************************************

      CALL get_nu_type(interaction,nu_type)

***************************************************************************

**** maximum differential cross section (without constant factor in front)
 
         IF (    nu_type .EQ. id_nu_e
     +      .OR. nu_type .EQ. id_nu_mu
     +      .OR. nu_type .EQ. id_nu_tau) THEN

           cv = 0.04
           ca = 1.27 / 2

         ELSE

           cv = 0.04
           ca = - 1.27 / 2

         ENDIF

         xsect_max = (cv + ca)**2 + (cv - ca)**2

***************************************************************************
*  Calculate differential cross section based on random cross section and
*  decide whether to keep cross section evaluated at recoil energy
***************************************************************************

***  Calculate maximum allowed electron kinetic energy

      T_p_max = 2.d0 * E_nu*E_nu / (2.d0 * E_nu + m_proton)

      done = .false.
         
      DOWHILE (.NOT. done)

         xsect_ran = ran_val(0.d0,1.d0)
         T_p = T_p_max * xsect_ran

* calculate test cross section

         xsect_diff = (cv + ca)**2  + (cv - ca)**2 * (E_nu - T_p)**2
     +                / E_nu**2  - (cv**2 - ca**2) * m_proton 
     +                * T_p / E_nu**2
         
         xsect_ratio = xsect_diff / xsect_max
     
         xsect_ran = ran_val(0.d0,1.d0)
         IF (xsect_ran .GT. xsect_ratio) THEN
            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO

*** Calculate proton kinetic energy

      E_part = T_p

C      xsect_diff = sig_e * xsect_diff
C      WRITE (*,*) E_part

***************************************************************************

      RETURN
      END
