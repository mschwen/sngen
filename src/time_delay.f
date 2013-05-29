************************************************************************
*  Filename: time_delay.f
************************************************************************
*
*  Modify time spectrum to include energy-dependent time delay effects for
*  massive neutrinos.  Even if all neutrinos are emitted at the same time
*  the dispersion experienced by massive neutrinos of different energies
*  is expected to stretch the burst in time.
*     
*  Inputs:
*       o  interaction ID (not array)
*	o  neutrino energy (MeV)
*	o  undelayed time (ms)
*  Ouputs:
*	o  delayed time for neutrino of specific mass and energy (ms)
*
************************************************************************

      REAL*8 FUNCTION time_delay(interaction,E_nu,time_undelay)

      IMPLICIT none

      INCLUDE 'physics_constants.inc'
      INCLUDE 'convert_units.inc'
      INCLUDE 'id_nu.inc'

C  Common Block Declarations  

      COMMON           /distance/
     +                 distance
      REAL*8           distance

      COMMON           /neutrino_masses/
     +                 nu_mass
      REAL*8           nu_mass(6)

C  Global Variable Declarations

      INTEGER          interaction
      DOUBLE PRECISION E_nu
      DOUBLE PRECISION time_undelay

C  Local Variable Declarations

      INTEGER          nu_type
      DOUBLE PRECISION m_nu, x,
     +                 delay_diff

**************************************************************************

***************************************************************************
*  Assign neutrino masses (eV/c^2)
***************************************************************************

* get neutrino type

      CALL get_nu_type(interaction,nu_type)

      IF (nu_type .EQ. id_nu_e) THEN
         m_nu = nu_mass(1)
      ELSEIF (nu_type .EQ. id_nu_ebar) THEN
         m_nu = nu_mass(2)
      ELSEIF (nu_type .EQ. id_nu_mu) THEN
         m_nu = nu_mass(3)
      ELSEIF (nu_type .EQ. id_nu_mubar) THEN
         m_nu = nu_mass(4)
      ELSEIF (nu_type .EQ. id_nu_tau) THEN
         m_nu = nu_mass(5)
      ELSEIF (nu_type .EQ. id_nu_taubar) THEN
         m_nu = nu_mass(6)
      ELSE
         WRITE (*,*) '*** time_delay.f: Something went wrong!'
         WRITE (*,*) '*** time_delay.f: --> did not assign mass'
      ENDIF

****************************************************************************

      IF (m_nu .GT. 0.d0) THEN

* convert neutrino mass units: eV/c^2 to MeV/c^2

      m_nu = m_nu * 1.d-6

***************************************************************************
*  Calculate the energy-dependent time delay experienced by a massive 
*  neutrino relative to a massless neutrino that travels at the speed of 
*  light in vacuum:
*    -> delay_diff = D(1/v - 1/c) = D/c[(1 - m^2c^4/E^2)^-1/2 - 1]
*    -> use binomial expansion: (1 + x)^-1/2   for -1 < x <= 1
*    -> (1 + x)^-1/2 ~ 1 - (1/2)x + ((1*3)/(2*4))x^2 - ((1*3*5)/(2*4*6))x^3
*       to the third order term where x = - (m^2c^4/E^2)
*    -> therefore (1 - (m^2c^4/E^2)^-1/2) ~ 1 + 1/2(m^2c^4/E^2)
*                                             - 3/8(m^2c^4/E^2)^2
*                                             + 15/48(m^2c^4/E^2)^3
*    -> hence delay_diff ~ D/c[1/2(m^2c^4/E^2) - 3/8(m^2c^4/E^2)^2
*                          + 15/48(m^2c^4/E^2)^3]
*
*    -> sadly, my thesis has an error: the last term in Eqns 3.6, 5.14
*       appears to have an extra factor of "2" in the denominator and is
*       thus incorrect :( (the "96" should read "48")
*
*  Units: [D] = kpc, [m_nu] = MeV, [E_nu] = MeV, [time_delay] = seconds
*
*  See Beacom & Vogel Phys Rev D, Vol 58, 093012 (1998) for simple formula:
*    -> delay_diff = 0.51463562514d0 * x * D  (to first order in the expansion)
*       where [D] = 10 kpc, [m_nu] = eV, [E_nu] = MeV, [time_delay] = seconds
***************************************************************************

         x = (m_nu * m_nu) / (E_nu * E_nu)

* energy-dependent time delay to third order (units of seconds)

         delay_diff = (distance * kpc2m / c_light)
     +               * (0.5d0 * x - 0.375d0 * x**2 + 0.3125d0 * x**3)

* convert time difference to milliseconds

         time_delay = time_undelay + (delay_diff / ms2s)
      ELSE
         time_delay = time_undelay
      ENDIF

      RETURN
      END
