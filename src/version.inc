***************************************************************************
*  Filename: version.inc
***************************************************************************
*
*  Version: 1.0 -- First "public" release (massless neutrinos only)
*  Version: 2.0 -- Include energy-dependent time delay for massive neutrinos
*           2.1 -- Multiply countrate data by number of supernovae
*               -- Use undelayed time for SN source number and spectrum 
*               -- Update physics constants to PDG 2002 values
*               -- Update targets (#d,#e) according to Robertson (ed) 2003
*               -- Update "O17_in_d2o" percentage (SNOMAN input)
*           2.2 -- Replace "ran1" routine in poidev.f with "ran_val"
*                  routine (which calls "ran3"); something similar is
*                  done in SNOMAN
*           2.3 -- Update physics constants to PDG 2003 values (ie., g_A)
*               -- Add code to calculate "central" energy values
*                  ( max[xsect(E) * fd_spectrum(E,time)] )
*               -- Consolodate "det_flux_es,cc,nc.f" into new "get_xsect"
*                  routine and remove det_flux_es,cc,nc.f routines
*           2.4 -- Edit group_sum2d_file to supress extra zero entries for
*                  unused delayed time bins
*               -- Update "en_nu_d_n_cc.f" to reflect changes made to similar
*                  SNOMAN routine "en_cc_neutron.for" (fixed cos_theta_n1 
*                  errors; round off solution in SNOMAN does not work in sngen,
*                  therefore keep "redo" flag to resample problem values)
*               -- Update physics constants to PDG 2004 values (ie., m_electron,*                  m_neutron,m_proton,m_nucleon,m_deuteron,m_n-m_p,G_fermi,
*                  cos_cabibbo (old one not correct 2002 value),sin^2(theta_W),
*                  mu_n,hbar_c2)
*           2.5 -- Include radiative corrections for nue,numu,nutau ES reactions
*                  (no published corrections for other ES reactions);
*                  -2.9% for nu_e ES counts, +1.8% for nu_mu, nu_tau ES counts;
*                  incorporated by spline fitting output values from QSNO
*               -- Update physics constant (M_A); affects nu_e_bar + p cross
*                  section (very minor)
*               -- Update total cross sections for all nu+d reactions 
*                  to NSAPMGK; count changes < 0.6% for each nu+d reaction
*               -- Update differential cross sections for nu+d CC reactions
*                  to NSAPMGK (leave nu+d NC with KN since new NC differential
*                  cross sections only go up to Enu=50 MeV; too low for SN);
*                  verified that XCSMFLAG (cross section smear flag) set to "F"
*           2.6 -- Revert to old "f_xs_diff" routine (2-point interpolation)
*                  and explicitly compile local copy (ie., include in Makefile);
*                  bug in RCON array proves fatal for 3-point interpolation for
*                  energies near cross section data step boundaries
*                  (eg., KN: 24.5, 54, 97 MeV); 2-point ok
*               -- Fix typo in Makefile (get_esdif_rad.o dependencies)
*               -- Revert to KN differential cross sections for sampling
*                  (give jagged final particle energy spectra, probably RCONS
*                  problem); continue to use NSAPMGK for count integration
*
****************************************************************************

      INTEGER*1        major_version, minor_version
      PARAMETER        (major_version = 2)
      PARAMETER        (minor_version = 6)

