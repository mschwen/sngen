***************************************************************************
* Filename: dir_nu_d_n_cc.f
***************************************************************************
*
*   Inputs:
*	o  neutrino energy (MeV)
*   Outputs:
*       o outgoing particles direction cosines
*
***************************************************************************

      SUBROUTINE dir_nu_d_n_cc(E_nu,part_dir)
   
      IMPLICIT none

      INCLUDE 'trig_constants.inc'
      INCLUDE 'physics_constants.inc'

C   Common Block Declarations

      COMMON /nu_d_positron_stuff/ E_e, sn_dir_cos, part_dir_e 
      REAL*8 E_e, sn_dir_cos(3), part_dir_e(3)

      COMMON /nu_d_E_neutrons/ E_n1, E_n2
      DOUBLE PRECISION E_n1, E_n2

      COMMON /nu_d_dir_n1/ cos_theta_n1
      REAL*8 cos_theta_n1

C   Global Variable Declarations

      REAL*8 E_nu, part_dir(3)

C   Local Variable Declarations

      REAL*8 dir_nu(3), dir_e(3)

      REAL*8 theta_n1

      DOUBLE PRECISION ran_val, azimuth_ran, azimuth

      REAL*8 p_e, p_n1_tot, p_n2_tot, p_n2(3)

      INTEGER ith_dir

      INTEGER dir_neutron_count
      DOUBLE PRECISION dir_n1(3)

      SAVE dir_neutron_count
      SAVE dir_n1

      DATA dir_neutron_count /0/
      DATA dir_n1 /3 * 0.d0/

***************************************************************************
*  Assign neutrino and positron directions
*  (information from COMMON block)
***************************************************************************

      dir_nu(1) = sn_dir_cos(1)
      dir_nu(2) = sn_dir_cos(2)
      dir_nu(3) = sn_dir_cos(3)

      dir_e(1) = part_dir_e(1)
      dir_e(2) = part_dir_e(2)
      dir_e(3) = part_dir_e(3)

****************************************************************************
****************************************************************************
*  Neutron #1 direction
*  [for nu_e_bar + d --> n + n + e+  CC reaction]
****************************************************************************

      IF (dir_neutron_count .EQ. 0) THEN      

         theta_n1 = dacos(cos_theta_n1)    ! determined from e+ kinematics
                                           ! (comes from COMMON block)

***  Choose the azimuthal scattering angle uniformly between 0 and 2*pi

         azimuth_ran = ran_val(0.d0,1.d0)
         azimuth = (2.d0 * pi) * azimuth_ran

***************************************************************************
*  Combine supernova neutrino angle with neutron #1 scattering angle
***************************************************************************

         CALL combine_angles(sn_dir_cos,theta_n1,azimuth,part_dir)

         DO 10 ith_dir = 1, 3
            dir_n1(ith_dir) = part_dir(ith_dir)
10       CONTINUE

***  Increment 'dir_neutron_count' for the next neutron in this event

         dir_neutron_count = 1

****************************************************************************
****************************************************************************
*  Neutron #2 direction
*  [for nu_e_bar + d --> n + n + e+  CC reaction]
****************************************************************************

      ELSEIF (dir_neutron_count .EQ. 1) THEN

* use conservation of momentum

***  positron momentum (lab frame)

        p_e = dsqrt(E_e**2 - m_electron**2)

***  neutron #1 and #2 momentum (lab frame)
***  (still something wrong here, I think)
 
        p_n1_tot = dsqrt(E_n1**2 - m_neutron**2)
c        p_n2_tot = dsqrt(E_n2**2 - m_neutron**2) ! this value is too small

        p_n2(1) = E_nu * dir_nu(1) - p_e * dir_e(1) 
     +           - p_n1_tot * dir_n1(1)
        p_n2(2) = E_nu * dir_nu(2) - p_e * dir_e(2) 
     +           - p_n1_tot * dir_n1(2)
        p_n2(3) = E_nu * dir_nu(3) - p_e * dir_e(3) 
     +           - p_n1_tot * dir_n1(3)

        p_n2_tot = dsqrt(p_n2(1)**2 + p_n2(2)**2 + p_n2(3)**2)

***  neutron #2 directional cosines

        part_dir(1) = p_n2(1) / p_n2_tot
        part_dir(2) = p_n2(2) / p_n2_tot
        part_dir(3) = p_n2(3) / p_n2_tot

c        WRITE (*,*) 'p_n2(1), p_n2(2), p_n2(3), p_n2_tot ', 
c     +               p_n2(1), p_n2(2), p_n2(3), p_n2_tot

***  Reset 'dir_neutron_count' back to '0' for the next neutrino event

        dir_neutron_count = 0

      ENDIF

      RETURN
      END
