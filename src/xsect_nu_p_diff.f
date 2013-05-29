****************************************************************************
*  Filename: xsect_nu_p_diff.f
****************************************************************************
*
*  This program defines the Llewellyn Smith differential quasi-elastic
*  cross section formalism,  d_(sigma) / d_(|q^2|), with the
*  modifications suggested by Vogel & Beacom.
*
*  See the following papers for more information:
*      o Llewellyn Smith, Physics Reports 3, 261-379 (1972)  [pgs 300-306]
*      o Vogel & Beacom hep-ph/9903554 April 1, 1999
*
*  For more recent constants than are to be found in Llewellyn Smith, see:
*      o Ahrens, et al. Physical Review  D35, 785 (1987)
*
*  Note: The results of this function have been compared with results 
*  supplied by J.F. Beacom (private communication); they are identical.
*
****************************************************************************

      FUNCTION xsect_nu_p_diff(E_nu,q2)

      IMPLICIT none

      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'

C   Global Variable Declarations

      DOUBLE PRECISION E_nu        ! input
      DOUBLE PRECISION q2          ! input  (integration variable)

C   Local Variable Declarations

      DOUBLE PRECISION xsect_nu_p_diff,   ! function
     +                 G_V_E, G_V_M,
     +                 F_V_1, xi_F_V_2, F_A,
     +                 sig_A, sig_B, sig_C,
     +                 s_minus_u,
     +                 dsig_dq2

      DOUBLE PRECISION m_e2, M, M2, m_n, m_p
      PARAMETER        (m_e2 = m_electron * m_electron)
      PARAMETER        (m_n = m_neutron)
      PARAMETER        (m_p = m_proton)
      PARAMETER        (M = m_nucleon)
      PARAMETER        (M2 = m_nucleon * m_nucleon)

****************************************************************************

***  form factors  [Llewellyn Smith equations: 3.16, 3.17, 3.24]

*     Sachs form factors:
      G_V_E = 1.d0 / ((1.d0 - (q2 / m_V2))**2)

      G_V_M = (1.d0 + xi) * G_V_E

*     Dirac electromagnetic isovector form factor:
      F_V_1 =  1.d0 / (1.d0 - (q2 / (4.d0 * M2))) 
     +        * (G_V_E - (q2 / (4.d0 * M2))
     +        * G_V_M)

*     Pauli electromagnetic isovector form factor:
      xi_F_V_2 = 1.d0 / (1.d0 - (q2 / (4.d0 * M2))) * (G_V_M - G_V_E)

*     Axial form factor:
      F_A = -g_A / (1.d0 - (q2 / m_A2))**2

***  cross section terms (in terms in form factors, etc.)
***  [Llewellyn Smith equations: 3.22]

      sig_A = (m_e2 - q2) / (4.d0 * M2)           ! factor out front
     +        * ((4.d0 - q2 / M2) * F_A**2      ! 1st term in brackets
     +        - (4.d0 + q2 / M2) * F_V_1**2     ! 2nd term in brackets
     +        - q2 / M2 * xi_F_V_2 * xi_F_V_2   ! 3rd term in brackets
     +         * (1.d0 + q2 / (4.d0 * M2))
     +        - 4.d0 * q2 * F_V_1 * xi_F_V_2 / M2
     +        - (m_e2 / M2) 
     +         * ((F_V_1 + xi_F_V_2)**2 
     +         + F_A * F_A))

      sig_B = -q2 / M2 * F_A * (F_V_1 + xi_F_V_2) 

      sig_C = 0.25d0 
     +       * (F_A * F_A
     +        + F_V_1 * F_V_1
     +        - q2 / (4.d0 * M2) 
     +        * xi_F_V_2 * xi_F_V_2)

* new Mandelstam variables  [Llewellyn Smith bottom of equation 3.18
*                           with modifications from Vogel & Beacom]
* -> s = (p_nu + p_p)**2 = m_p**2 + 2 * m_p * E_nu           (eqn A3)
* -> u = (p_nu - p_n)**2 = m_p**2 + m_e**2 - 2 * m_p * E_e   (eqn A5)

      s_minus_u =  4.d0 * m_p * E_nu + q2 - m_e2 - 2.d0 * M * delta_m

****************************************************************************
*  Differential cross section, d_(sigma) / d_(|q^2|)  
*  [Llewellyn Smith dsig_dq2 equation 3.18, with modifications from
*  Vogel & Beacom equation 16].
*
*  Multiply quasi-static cross section by 1.024 to account for inner
*  radiative corrections [see Vogel & Beacom equation 9].
*
*  dsig_dq2 units: cm^2
****************************************************************************

      dsig_dq2 = (M * G_fermi * cos_cabibbo)**2  
     +          * hbar_c2 / (8.d0 * pi * E_nu**2)
     +          * (sig_A
     +          + sig_B * s_minus_u / M2
     +          + (sig_C * s_minus_u**2
     +          - sig_C * 4.d0 * M2 * delta_m**2)
     +          / (M2 * M2))
     +          * 1.024d0

      xsect_nu_p_diff = dsig_dq2

      END
