***************************************************************************
* Filename: en_nu_p_cc.f
***************************************************************************
*
*   Inputs:
*	o  delayed neutrino interaction time (ms)
*	o  neutrino energy (MeV)
*   Outputs:
*	o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_p_cc(E_nu,id_part,E_part)
   
      IMPLICIT none

      INCLUDE 'id_part.inc'
      INCLUDE 'physics_constants.inc'

C   Global Variable Declarations

      INTEGER          id_part
      REAL*8           E_nu, E_part

C   Local Variable Declarations

      LOGICAL          done

      DOUBLE PRECISION a, b, c,
     +                 E_e_min, E_e_max,
     +                 q2_min, q2_max,
     +                 dsig_dcos_try(2), dsig_dcos_max

      DOUBLE PRECISION xsect_nu_p_diff,  ! functions
     +                 ran_val
      DOUBLE PRECISION q2_ran, dsig_dcos_ran,
     +                 xs_diff_ratio, xs_diff_ran

      DOUBLE PRECISION E_e, T_n

      DOUBLE PRECISION m_p, m_n, M, m_e2
      PARAMETER        (m_p = m_proton)
      PARAMETER        (m_n = m_neutron)
      PARAMETER        (M = m_nucleon)
      PARAMETER        (m_e2 = m_electron * m_electron)

***************************************************************************

* determine the range of E_e allowed from the quadratic equation
* (just depends on E_nu)

         a = (2.d0 * (E_nu + m_p))**2 - (2.d0 * E_nu)**2
         b = 2.d0 * (2.d0 * M * delta_m - 2.d0 * m_p * E_nu - m_e2) 
     +      * (2.d0 * (E_nu + m_p))
         c = (2.d0 * M * delta_m - 2.d0 * m_p * E_nu - m_e2)**2 
     +      + 4.d0 * m_e2 * E_nu**2

         E_e_min = (-b - dsqrt(b**2 - 4.d0*a*c))/(2.d0*a)
         E_e_max = (-b + dsqrt(b**2 - 4.d0*a*c))/(2.d0*a)

* for the chosen neutrino energy, E_nu, find maximum differential cross section
* value, dsig_dcos, to scale sample space
* -> either q2_min or q2_max results in highest differential cross section
*    (cross over point is E_nu ~20 MeV)

         q2_min = 2.d0*M*delta_m - 2.d0*m_p*(E_nu - E_e_min)
         q2_max = 2.d0*M*delta_m - 2.d0*m_p*(E_nu - E_e_max)

         dsig_dcos_try(1) = xsect_nu_p_diff(E_nu,q2_min)
         dsig_dcos_try(2) = xsect_nu_p_diff(E_nu,q2_max)

         IF (dsig_dcos_try(1) .GT. dsig_dcos_try(2)) THEN
            dsig_dcos_max = dsig_dcos_try(1)
         ELSE
            dsig_dcos_max = dsig_dcos_try(2)
         ENDIF

***************************************************************************
*  Sample differential cross section over the entire q^2 interval
***************************************************************************

      done = .false.

      DOWHILE (.NOT. done)

*** Calculate random d_(sigma) / d_(cos_theta)
*** (SN neutrino sampling has already taken care of the energy threshold)
      
         q2_ran = ran_val(q2_min,q2_max)

         dsig_dcos_ran = xsect_nu_p_diff(E_nu,q2_ran)
            
         xs_diff_ratio = dsig_dcos_ran / dsig_dcos_max

         xs_diff_ran = ran_val(0.d0,1.d0)
         IF (xs_diff_ran .GT. xs_diff_ratio) THEN
            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO

* convert chosen q2_ran into positron energy

      E_e = E_nu + q2_ran / (2.d0 * m_p)
     +     - 2.d0 * M * delta_m / (2.d0 * m_p)

*** positron (total) energy

      IF (id_part .EQ. id_positron) THEN

         E_part = E_e

*** neutron (kinetic) energy (using conservation of total energy)

      ELSEIF (id_part .EQ. id_neutron) THEN

         T_n = E_nu + m_p - (E_e + m_n)
         E_part = T_n

      ENDIF

      RETURN
      END
