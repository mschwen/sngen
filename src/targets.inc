*********************************************************************
*  Filename: targets.inc
**********************************************************************
*  The  media
**********************************************************************

      DOUBLE PRECISION N_h2o, N_d2o, N_acrylic, N_LAB, N_Pb

**** Number of light water molecules
**** 1.7 ktonne assuming 100% isotopic purity
      PARAMETER (N_h2o = 5.682749481d+31)

**** Number of heavy water molecules
      PARAMETER (N_d2o = 3.0377d+31)   ! Robertson (ed) May 2003

**** Number of acrylic molecules
**** 34 tonne acrylic vessel: CH2=C(CH3)CO2CH3  or C5H8O2 for short
*        o  33 tonne acrylic vessel (including thicker upper hemisphere
*           and rope groove panels at the equator)
*        o  1 tonne chimney visible by PMTs (2.5 tonnes for total length)
      PARAMETER (N_acrylic = 2.04515776d+29)

**** Number of LAB molecules
**** 860 tons of LAB: C12-C6H6
      PARAMETER (N_LAB = 2.10d30)

**** Number of Pb molecules
**** 80 tons of Pb208
      PARAMETER (N_Pb = 2.90d29)

*********************************************************************
*  The specific SNO targets
*********************************************************************

      DOUBLE PRECISION N_e_h2o, N_e_d2o
      DOUBLE PRECISION N_proton_h2o, N_proton_av, N_proton_d2o
      DOUBLE PRECISION N_deuteron
      DOUBLE PRECISION N_12C_AV
      DOUBLE PRECISION N_16O_h2o, N_16O_d2o, N_16O_av 
      DOUBLE PRECISION N_18O_h2o, N_18O_d2o
      DOUBLE PRECISION N_e_LAB, N_proton_LAB, N_12C_LAB

      PARAMETER (N_e_h2o = (8.d0 + 2.d0) * N_h2o)   ! from all molecules
      PARAMETER (N_e_d2o = 30.4046d+31)             ! Robertson (ed) May 2003
*      PARAMETER (N_e_d2o = (8.d0 + 2.d0) * N_d2o)   ! from all molecules
      PARAMETER (N_proton_h2o = 2.d0 * N_h2o)       ! only from H
      PARAMETER (N_proton_d2o = 2.d0 * (8.24d-04 * N_d2o)) ! h2o in d2o      
      PARAMETER (N_proton_av = 8.d0 * N_acrylic)    ! only from H
      PARAMETER (N_deuteron = 2.d0 * N_d2o)
      
      PARAMETER (N_12C_AV = 5.d0 * N_acrylic)
      PARAMETER (N_16O_h2o = N_h2o)
      PARAMETER (N_16O_d2o = N_d2o)
      PARAMETER (N_16O_av = 2.d0 * N_acrylic)
      PARAMETER (N_18O_h2o = 2.04d-03 * N_h2o)
      PARAMETER (N_18O_d2o = 2.04d-03 * N_d2o)

      PARAMETER (N_e_LAB = (138.d0) * N_LAB) ! from all molecules
      PARAMETER (N_proton_LAB = (30.d0) * N_LAB)       ! only from H
      PARAMETER (N_12C_LAB = 18.d0 * N_LAB)
