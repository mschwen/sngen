***************************************************************************
* Filename: dir_nu_p_cc.f
***************************************************************************
*
*   Inputs:
*	o  total neutrino energy (MeV)
*	o  total positron energy  (MeV)
*	o  supernova direction cosines
*   Outputs:
*	o  sampled final state particle directional cosines
*
***************************************************************************

      SUBROUTINE dir_nu_p_cc(E_nu,id_part,E_part,sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'id_part.inc'
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'

C   Global Variable Declarations

      INTEGER          id_part                 ! input
      REAL*8           E_nu, E_part            ! input
      REAL*8           sn_dir_cos(*)           ! input
      REAL*8           part_dir(3)             ! output

C   Local Variable Declarations

      DOUBLE PRECISION p_e, E_n, p_n

      DOUBLE PRECISION ran_val
      DOUBLE PRECISION cos_theta, theta, azimuth

      DOUBLE PRECISION m_electron2, m_proton2, m_neutron2
      PARAMETER        (m_electron2 = m_electron * m_electron)
      PARAMETER        (m_proton2 = m_proton * m_proton)
      PARAMETER        (m_neutron2 = m_neutron * m_neutron)

***************************************************************************

      IF (id_part .EQ. id_positron) THEN

* SNOMAN formula agrees with Beacom (private communication)

         p_e = dsqrt(E_part**2 - m_electron2)
         cos_theta = (0.5d0 * (m_neutron2 - m_proton2 - m_electron2)
     +              + E_part * m_proton + E_nu * (E_part - m_proton))
     +              / (E_nu * p_e)

      ELSEIF (id_part .EQ. id_neutron) THEN

         E_n = E_part + m_neutron
         p_n = dsqrt(E_n**2 - m_neutron**2)
         cos_theta = (m_proton * (E_n - E_nu) + E_nu * E_n
     +              - 0.5d0 * (m_neutron2 + m_proton2 - m_electron2))
     +              / (E_nu * p_n)

      ENDIF

* if 'cos_theta' is outside the bounds, set it to +/- 1 (depending on what
* the sign of the original value was)

      IF (abs(cos_theta) .GT. 1.d0) THEN
         cos_theta = sign(1.d0,cos_theta)
      ENDIF

      theta = dacos(cos_theta)

* Choose the azimuthal scattering angle uniformly between 0 and 2*pi

      azimuth = (2.d0 * pi) * ran_val(0.d0,1.d0)

***************************************************************************
*  Combine supernova angle with the final state particle scattering angle
***************************************************************************

      CALL combine_angles(sn_dir_cos,theta,azimuth,part_dir)

***************************************************************************

      RETURN
      END
