***************************************************************************
*  Filename: threshold.inc
***************************************************************************

      DOUBLE PRECISION thres_nu_e_bar_p_cc,
     +                 thres_nu_e_d_cc, thres_nu_e_bar_d_cc,
     +                 thres_nu_d_nc, thres_nu_16O_nc, 
     +                 thres_nu_e_O_cc, thres_nu_e_bar_O_cc,
     +                 thres_nu_12C_nc,
     +                 thres_nu_e_16O_cc, thres_nu_e_bar_16O_cc, 
     +                 thres_nu_e_12C_cc, thres_nu_e_bar_12C_cc,
     +                 thres_nu_es_e, thres_nu_es_p,
     +                 thres_nu_208Pb_cc_1n, thres_nu_208Pb_cc_2n,
     +                 thres_nu_208Pb_nc_1n, thres_nub_208Pb_nc_1n,
     +                 thres_nu_208Pb_nc_2n, thres_nub_208Pb_nc_2n

      PARAMETER (thres_nu_es_e = 0.d0)
      PARAMETER (thres_nu_es_p = 0.d0)
      PARAMETER (thres_nu_e_d_cc = 1.44278739699d0)     ! using NSGK formula
                                                        ! (PDG 2000 constants)
      PARAMETER (thres_nu_d_nc = 2.224566d0)            ! m_n + m_p - m_d
                                                        ! (PDG 2000 constants)
      PARAMETER (thres_nu_e_12C_cc = 17.30d0)
      PARAMETER (thres_nu_e_bar_12C_cc = 13.40d0)
      PARAMETER (thres_nu_e_O_cc = 1.66d0)              ! set to 18O thres
      PARAMETER (thres_nu_e_16O_cc = 15.4d0)
      PARAMETER (thres_nu_e_bar_16O_cc = 10.4d0)
      PARAMETER (thres_nu_e_bar_O_cc = 11.40d0)
      PARAMETER (thres_nu_12C_nc = 15.11d0)
      PARAMETER (thres_nu_16O_nc = 15.00d0)
      PARAMETER (thres_nu_e_bar_p_cc = 1.80606579874d0) ! using Beacom formula
                                                        ! (PDG 2000 constants)
      PARAMETER (thres_nu_e_bar_d_cc = 4.03322402394d0) ! using NSGK formula
                                                        ! (PDG 2000 constants)
      PARAMETER (thres_nu_208Pb_cc_1n = 13.0)  ! Estimated from PhysRevD67,013005
      PARAMETER (thres_nu_208Pb_cc_2n = 25.0)  
      PARAMETER (thres_nu_208Pb_nc_1n = 11.0)  
      PARAMETER (thres_nub_208Pb_nc_1n = 11.0)  
      PARAMETER (thres_nu_208Pb_nc_2n = 26.0)
      PARAMETER (thres_nub_208Pb_nc_2n = 26.0)

