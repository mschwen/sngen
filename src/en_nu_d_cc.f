***************************************************************************
* Filename: en_nu_d_cc.f
***************************************************************************
*
*  This routine samples the cross section in order to obtain acceptable
*  electron/positron energies for the the reactions:
*       *  nu_e + d --> p + p + e-  CC (D2O)
*       *  nu_e_bar + d --> n + n + e+  CC (D2O)
*
*  Similar to SNOMAN routine "en_cc.for"
*
*  Inputs:
*	o  interaction id number (not array)
*	o  delayed neutrino interaction time (ms)
*	o  neutrino energy (MeV)
*  Outputs:
*	o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_d_cc(interaction,E_nu,E_part)
   
      IMPLICIT none

      INCLUDE 'physics_constants.inc'

C   Common Block Declarations

      INTEGER int_type ! needed by 'f_xs_diff.f'
      REAL Enu         ! needs to be REAL for 'f_xs_diff.f' not REAL*8
      COMMON /nu_d_en_stuff/ Enu, int_type
      REAL*8 cos_theta_e ! needed by 'dir_nu_d_cc.f' and 'dir_nu_d_n_cc.f'
      COMMON /nu_d_dir_e/ cos_theta_e

C   Global Variable Declarations

      INTEGER interaction
      REAL*8 E_nu, E_part

C   Local Variable Declarations

      REAL*8 E_thres
      REAL*8 T_e_max

      EXTERNAL xsect_nu_d_diff
      REAL p(2), x_i(2,2), iter, fret, xsect_nu_d_diff ! REAL for 'powell.f'

      REAL*8 xsect_max, xsect_diff, xsect_ratio

      LOGICAL done


      DOUBLE PRECISION ran_val
      DOUBLE PRECISION energy_ran
      DOUBLE PRECISION cos_ran
      REAL T_e_ran, E_e_ran, cos_theta_ran  ! needs to be REAL for 'f_xs_diff.f'
      DOUBLE PRECISION xsect_ran

      REAL*8 f_xs_diff

***************************************************************************
*  Get reaction interaction id
***************************************************************************

      CALL get_int_type(interaction,int_type)
      Enu = E_nu

***************************************************************************
*  Calculate maximum electron (or positron) kinetic energy, 'T_e_max'
***************************************************************************

      CALL get_threshold(interaction,E_thres)
      T_e_max = E_nu - E_thres

***************************************************************************
*  Find maximum of differential cross section for this neutrino energy,
*  'E_nu', using Powell's method
***************************************************************************

      p(1) = -0.999          ! cos_theta starting values
      p(2) = 0.002 * E_nu    ! E_nu starting value
      x_i(1,1) = 1.0
      x_i(1,2) = 0.0
      x_i(2,1) = 0.0
      x_i(2,2) = 1.0

      CALL powell(p,x_i,2,2,0.001,iter,fret,xsect_nu_d_diff)
      
      xsect_max = -fret

***************************************************************************
*  Sample electron energy and electron scattering angle, 'theta'
***************************************************************************

      done = .false.
         
      DOWHILE (.NOT. done)

* get random electron/positron kinetic energy value

         energy_ran = ran_val(0.d0,1.d0)
         T_e_ran = T_e_max * energy_ran ! KN cross sections require kinetic E
c         E_e_ran = T_e_ran + m_electron ! NSAPMGK cross sections require total E

* get random 'cos_theta' value (no constraint from above kinetic energy)

         cos_ran = ran_val(-1.d0,1.d0)
         cos_theta_ran = cos_ran      ! convert to REAL value for SNOMAN

* call SNOMAN 'nu + d' differential cross section routine
* (now takes total energy)

         xsect_diff = f_xs_diff(int_type,Enu,cos_theta_ran,T_e_ran) ! KN
c         xsect_diff = f_xs_diff(int_type,Enu,cos_theta_ran,E_e_ran) ! NSAPMGK

         xsect_ratio = xsect_diff / xsect_max

* sample energy and theta at the same time by using the doubly
* differential cross section
         xsect_ran = ran_val(0.d0,1.d0)
         IF (xsect_ran .GT. xsect_ratio) THEN
            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO

*************************************************************************** 
*  Now that we've finished sampling, calculate the total energy for the
*  final state electron (or positron) and the electron (or positron)
*  scattering angle, theta
*
*  This information will be used for both the nu_e + d CC reaction and the
*  nu_e_bar + d CC reaction
***************************************************************************

      E_part = T_e_ran + m_electron

      cos_theta_e = cos_ran

      RETURN
      END
