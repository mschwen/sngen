***************************************************************************
* Filename: en_nu_d_nc.f
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

      SUBROUTINE en_nu_d_nc(interaction,E_nu,E_part)
   
      IMPLICIT none

      INCLUDE 'physics_constants.inc'

C   Common Block Declarations

      INTEGER int_type ! needed by 'f_xs_diff.f'
      REAL Enu         ! needs to be REAL for 'f_xs_diff.f' not REAL*8
      COMMON /nu_d_en_stuff/ Enu, int_type
      REAL*8 costheta  ! needed by 'dir_nu_d_nc.f'
      COMMON /nu_d_dir_stuff/ costheta

C   Global Variable Declarations

      INTEGER interaction
      REAL*8 E_nu, cos_theta, E_part

C   Local Variable Declarations

      REAL*8 T_n_max, E_n_max
      REAL*8 s
      REAL*8 beta, gamma
      REAL*8 p_n_cms, E_n_cms

      EXTERNAL xsect_nu_d_diff
      REAL p(2), x_i(2,2), iter, fret, xsect_nu_d_diff

      REAL*8 xsect_max, xsect_diff, xsect_ratio

      LOGICAL done

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION energy_ran
      DOUBLE PRECISION cos_ran
      REAL T_n_ran, cos_theta_ran     ! needs to be REAL for 'f_xs_diff.f'
      DOUBLE PRECISION xsect_ran

      REAL*8 f_xs_diff

***************************************************************************
*  Set up common block variables for 'xsect_nu_d_diff.f' since we can only
*  pass one argument to the function because of how it is used in the
*  'powell.for' maximization search routine
***************************************************************************

      CALL get_int_type(interaction,int_type)
      Enu = E_nu

***************************************************************************
*  Calculate maximum neutron kinetic energy, 'T_n_max'
***************************************************************************

***  total energy of neutrino-deuteron system
      s = m_deuteron**2 + 2.d0 * E_nu * m_deuteron

***  momentum of neutron in CMS frame

      p_n_cms = dsqrt((s - (m_proton + m_neutron)**2)
     +         * (s - delta_m**2)) / (2.d0 * dsqrt(s))

***  total energy of neutron in CMS frame

      E_n_cms = dsqrt(p_n_cms**2 + m_neutron**2)

*** velocity of CMS in the lab frame: PDG pg. 250 Eqn (34.4)
*   (p_nu = E_nu)

      beta = E_nu / (E_nu + m_deuteron)
      gamma = 1.d0 / dsqrt(1.d0 - beta**2)

*** maximum neutron energy in the lab

      E_n_max = gamma * (beta * p_n_cms + E_n_cms)

***  kinetic energy of neutron (with correction for nucleon motion)

      T_n_max = E_n_max - m_neutron
      T_n_max = T_n_max * 1.05d0     ! 5% correction for nucleon motion
 
***************************************************************************
* Determine maximum of differential cross section for this neutrino
* energy using Powell's method
***************************************************************************

      p(1) = -0.999          ! cos_theta starting value
      p(2) = 0.002 * E_nu    ! E_nu starting value
      x_i(1,1) = 1.0
      x_i(1,2) = 0.
      x_i(2,1) = 0.
      x_i(2,2) = 1.0

      CALL powell(p,x_i,2,2,0.001,iter,fret,xsect_nu_d_diff)
      
      xsect_max = -fret

***************************************************************************
*  Sample neutron kinetic energy and neutron scattering angle, 'theta'
***************************************************************************

      done = .false.
         
      DOWHILE (.NOT. done)
         energy_ran = ran_val(0.d0,1.d0)
         T_n_ran = T_n_max * energy_ran

         cos_ran = ran_val(-1.d0,1.d0)
         cos_theta_ran = cos_ran

*** call SNOMAN 'nu + d' differential cross section routine
*** (both KN and NSAPMGK expect neutron kinetic energy, so no need to change)
         xsect_diff = f_xs_diff(int_type,Enu,cos_theta_ran,T_n_ran)

c         WRITE (*,*) 'E_nu, cos_theta_ran, E_e_ran, xsect_diff =', 
c     +                Enu,cos_theta_ran,E_e_ran,xsect_diff

         xsect_ratio = xsect_diff / xsect_max
     
         xsect_ran = ran_val(0.d0,1.d0)
         IF (xsect_ran .GT. xsect_ratio) THEN
            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO
 
***  Return neutron kinetic energy 

      E_part = T_n_ran

      cos_theta = cos_theta_ran
      costheta = cos_theta

      RETURN
      END
