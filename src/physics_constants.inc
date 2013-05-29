****************************************************************************
*  Filename: physics_constants.inc
****************************************************************************
*
*  Consolidation of physics constants used in the SNO supernova
*  generator, 'sngen'
*
****************************************************************************

      DOUBLE PRECISION c_light,
     +                 m_electron,
     +                 m_proton, m_neutron, delta_m,
     +                 m_nucleon,
     +                 m_deuteron,
     +                 m_V2, m_A2,
     +                 G_fermi,
     +                 cos_cabibbo, sin2theta,
     +                 xi, g_A,
     +                 hbar_c2,
     +                 alpha
c     +                 E_swap_nu,
c     +                 E_swap_nu_bar

**** [speed of light in vacuum: Particle Data Group 2004  (m/s)]
      PARAMETER (c_light = 2.99792458d+08)

**** [electron mass: Particle Data Group 2004 Physical Constants  (MeV/c^2)]
      PARAMETER (m_electron = 0.510998918d0)

**** [proton mass: Particle Data Group 2004 Physical Constants  (MeV/c^2)]
      PARAMETER (m_proton = 938.272029d0)

**** [neutron mass: Particle Data Group 2004 Baryon Properties  (MeV/c^2)]
      PARAMETER (m_neutron = 939.56536d0)

**** [average nucleon mass using the above values  (MeV/c^2)]
      PARAMETER (m_nucleon = 938.918695d0)

**** [deuteron mass: Particle Data Group 2004  (MeV/c^2)]
      PARAMETER (m_deuteron = 1875.61282d0)

**** [vector dipole mass squared:
*     Ahrens et al. Phys Rev D 35 785 (1987) m_V = 840 MeV/c^2  (MeV/c^2)^2]
      PARAMETER (m_V2 = 7.056d05)

**** [axial-vector dipole mass squared: 
*     Bernard et al. J. Phys G 28 R1 (2002) m_A = 1026 MeV/c^2  (MeV/c^2)^2]
      PARAMETER (m_A2 = 1.052676d06)

**** [nucleon mass difference (m_n-m_p): Particle Data Group 2004 
*     Baryon Properties  (MeV/c^2)]
      PARAMETER (delta_m = 1.2933317d0)

**** [Fermi coupling constant: Particle Data Group 2004  (MeV)^-2]
      PARAMETER (G_fermi = 1.16637d-11)

**** [cos(Cabibbo angle): Particle Data Group 2004 |V_ud|=0.9738]
      PARAMETER (cos_cabibbo = 0.9738d0)

**** [weak mixing angle squared, sin^2(theta_W): Particle Data Group 2004]
      PARAMETER (sin2theta = 0.23120d0)  ! M_Z (MS-bar) value a la Bahcall

**** [anomolous magnetic moment, xi = (mu_p - 1) - mu_n: 
*     Particle Data Group 2004: (mu_p - 1) = 1.792847337; mu_n = -1.9130427]
      PARAMETER (xi = 3.705890037d0)

**** [weak axial coupling constant: PDG 2004 (neutron decay parameters)]
      PARAMETER (g_A = 1.2695d0)

**** [(h_bar * c)^2: Particle Data Group 2004  (MeV^2 * cm^2)]
      PARAMETER (hbar_c2 = 3.89379323d-22)

**** [alpha = 1/137.xxx: Particle Data Group 2004
      PARAMETER (alpha = 7.297352568d-3)

**** [E_swap_nu: Neutrino energy above which flavour swapping occurs: JCPA12(2007)010]
c      PARAMETER (E_swap_nu = 7)
c      PARAMETER (E_swap_nu = 0)

**** (E_swap_nu_bar: JCPA12(2007)010]
c      PARAMETER (E_swap_nu_bar = 4)
