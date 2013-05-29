***********************************************************************
*  Filename: read_user.f
***********************************************************************
*
* Loads all of the constants from "sn_constants.dat" into the 
* appropriate COMMON blocks
*
**********************************************************************

      SUBROUTINE read_user()

      IMPLICIT none

      INCLUDE 'file_flux_size.inc'
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      REAL*8           distance
      COMMON           /distance/ 
     +                 distance

      REAL*8           nu_mass(6)
      COMMON           /neutrino_masses/
     +                 nu_mass

      DOUBLE PRECISION declination, right_ascension
      DOUBLE PRECISION dsec
      INTEGER          imin, ihr, iday, imonth, iyear
      COMMON           /direction/ 
     +                 declination, right_ascension,
     +                 dsec, imin, ihr, iday, imonth, iyear

      INTEGER          sn_model_switch
      COMMON           /sn_model/
     +                 sn_model_switch
   
      INTEGER          swap_flag
      COMMON           swap_flag

* Burrows
      REAL*8           lts(3,0:9), ltau(3,6), lco(3,11), lsig(3,4),
     +                 lw(3), ltheta(3), ltn(3), ldelt(3),
     +                 ia(3), ito(3), ib(3),
     +                 bc(3,6), bta(3),
     +                 btb(3), bsiga(3), bsigb(3), bn(3), bdt(3)
      COMMON           /Burrows_lum_stuff/
     +                 lts, ltau, lco, lsig, lw, ltheta, ltn, ldelt
      COMMON           /Burrows_eta_stuff/
     +                 ia, ito, ib
      COMMON           /Burrows_temp_stuff/
     +                 bc, bta, btb, bsiga, bsigb, bn, bdt

* Beacom&Vogel
      REAL*8           E_bind, tau1, tau2,
     +                 temp_nu_e, temp_nu_eb, temp_nu_mu
      COMMON           /BeacomVogel_lum_stuff/
     +                 E_bind, tau1, tau2
      COMMON           /BeacomVogel_temp_stuff/
     +                 temp_nu_e, temp_nu_eb, temp_nu_mu

      INTEGER          use_author_params
      COMMON           /sn_model_params/
     +                 use_author_params

      INTEGER          model_stats
      COMMON           /sn_model_stats/
     +                 model_stats

      INTEGER          n_supernova
      COMMON           /number_supernova/
     +                 n_supernova

      INTEGER          gen_flag
      CHARACTER*40     counts_file
      COMMON           /generator_flag/
     +                 gen_flag, counts_file

c      CHARACTER*80     lumen_file, fermi_file
c      COMMON           /sn_spectra_files/
c     +                 lumen_file, fermi_file

      INTEGER          lumen_fit_file_flag
      CHARACTER*40     lumen_nu_e_file, lumen_nu_eb_file, 
     +                 lumen_nu_mu_file
      COMMON           /sn_lumen_fit_file/
     +                 lumen_fit_file_flag, 
     +                 lumen_nu_e_file, lumen_nu_eb_file,
     +                 lumen_nu_mu_file

      INTEGER          lum_cutoff_flag
      REAL*8           lum_cutoff_time
      COMMON           /lum_cutoff_stuff/
     +                 lum_cutoff_time,
     +                 lum_cutoff_flag

      INTEGER          FD_shape_file_flag
      REAL*8           FD_time
      CHARACTER*40     FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu
      COMMON           /FD_shape/
     +                 FD_time, FD_shape_file_flag,
     +                 FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu

      INTEGER          central_flag
      CHARACTER*40     central_file
      COMMON           /central_stuff/
     +                 central_flag, central_file

      REAL*8           emin, emax, estep
      INTEGER          dum_emin, dum_emax, dum_estep
      COMMON           /energy_limits/
     +                 emin, emax, estep,
     +                 dum_emin, dum_emax, dum_estep

      REAL*8           time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2
      INTEGER          dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2
      COMMON           /time_limits/
     +                 time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2,
     +                 dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2

***  Calculated here, not read in from user file
      INTEGER          energy_entries
      COMMON           /sn_energy_entries/
     +                 energy_entries

***  Calculated here, not read in from user file
      INTEGER          time_entries1, time_entries2, time_entries 
      COMMON           /sn_time_entries/
     +                 time_entries1, time_entries2, time_entries

***  Calculated here, not read in from user file
      INTEGER          flux_entries1, flux_entries2, flux_entries
      COMMON           /sn_flux_entries/
     +                 flux_entries1, flux_entries2, flux_entries

      INTEGER          xsection_file_flag
      CHARACTER*40     xsection_file
      COMMON           /sn_xsection_file/
     +                 xsection_file_flag, xsection_file

      REAL*8           T_e_min
      COMMON           /det_threshold/
     +                 T_e_min

      CHARACTER*40     all_crfile
      CHARACTER*40     es_crfile, cc_crfile, nc_crfile
      CHARACTER*40     nu_e_crfile, nu_eb_crfile, nu_mu_crfile
      CHARACTER*40     d2o_crfile, h2o_crfile, av_crfile
      COMMON           /nu_countrate_file/
     +                 all_crfile,
     +                 es_crfile, cc_crfile, nc_crfile,
     +                 nu_e_crfile, nu_eb_crfile, nu_mu_crfile,
     +                 d2o_crfile, h2o_crfile, av_crfile

      CHARACTER*40     all_cfile
      CHARACTER*40     es_cfile, cc_cfile, nc_cfile
      CHARACTER*40     nu_e_cfile, nu_eb_cfile, nu_mu_cfile
      CHARACTER*40     d2o_cfile, h2o_cfile, av_cfile
      COMMON           /nu_count_file/
     +                 all_cfile,
     +                 es_cfile, cc_cfile, nc_cfile,
     +                 nu_e_cfile, nu_eb_cfile, nu_mu_cfile,
     +                 d2o_cfile, h2o_cfile, av_cfile

      CHARACTER*40     all_cumcfile
      CHARACTER*40     es_cumcfile, cc_cumcfile, nc_cumcfile
      CHARACTER*40     nu_e_cumcfile, nu_eb_cumcfile, nu_mu_cumcfile
      CHARACTER*40     d2o_cumcfile, h2o_cumcfile, av_cumcfile
      COMMON           /nu_cumulativecount_file/
     +                 all_cumcfile,
     +                 es_cumcfile, cc_cumcfile, nc_cumcfile,
     +                 nu_e_cumcfile, nu_eb_cumcfile, nu_mu_cumcfile,
     +                 d2o_cumcfile, h2o_cumcfile, av_cumcfile

      INTEGER pois_flag
      COMMON /sn_pois_flag/
     + pois_flag

      REAL*8 trigger_gate
      INTEGER time_overlap_flag
      CHARACTER*40 time_overlap_file
      COMMON /sn_time_overlap/
     + trigger_gate, time_overlap_flag, time_overlap_file

      CHARACTER*40 nu_sample_file
      COMMON /nu_sample_file/
     + nu_sample_file

      CHARACTER*40 part_sample_file
      COMMON /part_sample_file/
     + part_sample_file

      INTEGER efficiency_flag
      COMMON /sample_efficiency/
     + efficiency_flag

      INTEGER mcpl_flag
      CHARACTER*40 mcpl_file, mcpl_default
      COMMON /mcpl_stuff/
     + mcpl_flag, mcpl_file, mcpl_default

      INTEGER cmd_flag
      CHARACTER*40 cmd_default, zdab_default, hbook_default
      INTEGER zdab_flag
      COMMON /cmd_stuff/
     + cmd_flag, cmd_default, zdab_flag, zdab_default, hbook_default

      REAL*8 nhit_thres, comp_thres,
     + pmt_efficiency, h2o_in_d2o, O17_in_d2o
      COMMON /snoman_parameters/
     + nhit_thres, comp_thres, pmt_efficiency, h2o_in_d2o, O17_in_d2o

C  Local Variable Declarations

      INTEGER          nu
      REAL*8           sn_time
      CHARACTER*80     description

****************************************************************************

****************************************************************************
*   Open file containing SN constants
****************************************************************************

      OPEN (unit=10, file='dat/sn_user.dat', status='old')

****************************************************************************
*   Read in constants from "sn_user.dat" file
****************************************************************************

      WRITE (*,*) 'Loading SN generator user input parameters...'
      WRITE (*,*) 

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

C  Distance to SN (in kpc)

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1000) distance

      sn_time = ((distance * kpc2m) / c_light) * sec2year

      WRITE (*,*) 'A supernova exploded at',distance,' kpc',sn_time,
     +' years ago!'
      WRITE (*,*)

C  Neutrino Masses (in eV/c^2)

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1000) nu_mass(1)
      READ (10,9000) description
      READ (10,1000) nu_mass(2)
      READ (10,9000) description
      READ (10,1000) nu_mass(3)
      READ (10,9000) description
      READ (10,1000) nu_mass(4)
      READ (10,9000) description
      READ (10,1000) nu_mass(5)
      READ (10,9000) description
      READ (10,1000) nu_mass(6)

C  SN coordinates in "celestial" reference frame (degrees, degrees, hours)

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1000) declination  
      READ (10,9000) description
      READ (10,1000) right_ascension
      READ (10,9000) description
      READ (10,1099) iyear, imonth, iday, ihr, imin, dsec

C  Supernova Model

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,9000) description
      READ (10,1001) sn_model_switch 
      READ (10,9000) description
      READ (10,1001) use_author_params
      READ (10,9000) description
      READ (10,1001) model_stats

*****************************************************************************
* Burrows et al. (1992)
*****************************************************************************
      IF (sn_model_switch .EQ. 1) THEN
* energy (MeV)
         emin  = 1.00000d-01
         emax  = 100.0000d+00
         estep = 1.00000d-01
         energy_entries = nint((emax - emin) / estep + 1)
* time (ms)
         tstep1      = 1.00000d+00    ! use 1ms steps for "burst"
c         tstep1      = 1.00000d+01    ! use 10ms steps for "burst"
         time_start1 = -1.00000d+02
         time_stop1  = 1.00000d+03    ! out to 1 second for "burst"
c         tstep2      = 1.00000d+01    ! use 10ms steps for "cooling"
         tstep2      = 2.00000d+01    ! use 20ms steps for "cooling" Mac OS X
c         tstep2      = 1.00000d+02    ! use 100ms steps for "cooling"
         time_start2 = time_stop1 + tstep2    ! tstep2 off from time_stop1
         time_stop2  = 5.00000d+04    ! out to 50 seconds for "cooling"

         time_entries1 = nint((time_stop1 - time_start1) / tstep1 + 1)
         time_entries2 = nint((time_stop2 - time_start2) / tstep2 + 1)
         time_entries = time_entries1 + time_entries2

* read Burrows model parameters
         OPEN (unit = 13, file = 'dat/model/Burrows/Burrows.dat',
     +         status = 'old')
         READ (13,*) description
         READ (13,*) description
         READ (13,*) description
         READ (13,*) description
         READ (13,*) description
         READ (13,*) description
         READ (13,*) description
         DO nu = 1, 3       ! 3 general neutrino types
            READ (13,*) description
            READ (13,*) description,LTS(nu,1),LTAU(nu,1)
            READ (13,*) description,LTS(nu,2),LTAU(nu,2)
            READ (13,*) description,LTS(nu,3),LCO(nu,3),LTS(nu,0),
     +                  LSIG(nu,1)
            READ (13,*) description,LTS(nu,4),LTAU(nu,3)
            READ (13,*) description,LTS(nu,5),LCO(nu,5),LTAU(nu,4),
     +                  LCO(nu,6), LTAU(nu,5),LW(nu),LTHETA(nu)
            READ (13,*) description,LTS(nu,7),LCO(nu,7),LTS(nu,6)
            READ (13,*) description,LSIG(nu,3)
            READ (13,*) description,LTAU(nu,6)
            READ (13,*) description,LTS(nu,9),LCO(nu,10),LTN(nu),
     +                  LDELT(nu)
            READ (13,*) description,LSIG(nu,4)
            READ (13,*) description,IA(nu),ITO(nu),IB(nu)
            READ (13,*) description,BC(nu,1),BTA(nu)
            READ (13,*) description,BC(nu,2),BC(nu,3)
            READ (13,*) description,BC(nu,4),BTB(nu),BSIGA(nu),
     +                  BSIGB(nu)
            READ (13,*) description,BC(nu,6),BN(nu),BDT(nu)

* determine needed constants to make lumonsity curve continous
 
            LW(nu) = (2.d0 * pi) / LW(nu)
            LTHETA(nu) = LTHETA(nu)*deg2rad - LW(nu)*LTS(nu,4)
            LCO(nu,6) = LCO(nu,6) / dexp(LTS(nu,4) / LTAU(nu,5))
            IF (LTAU(nu,2) .NE. 0.d0) THEN
               LCO(nu,2) = LCO(nu,3)
     +                   * dexp((LTS(nu,2) - LTS(nu,0))**2
     +                   * (-0.5d0) / LSIG(nu,1)**2
     +                   - LTS(nu,2) / LTAU(nu,2))
            ELSE
               LCO(nu,2)  = 0.d0
               LTAU(nu,2) = 1.d0
            END IF
            LCO(nu,1) = LCO(nu,2)
     +                * dexp(LTS(nu,1) / LTAU(nu,2) 
     +                - LTS(nu,1) / LTAU(nu,1))
            LCO(nu,4) = LCO(nu,3)
     +                * dexp((LTS(nu,3) - LTS(nu,0))**2
     +                * (-0.5d0) / LSIG(nu,1)**2
     +                - LTS(nu,3) / LTAU(nu,3))
            IF (nu .EQ. 1)
     +       LCO(nu,5) = (LCO(nu,4)
     +                 * dexp(LTS(nu,4) / LTAU(nu,3))
     +                 - LCO(nu,6)
     +                 * dexp(LTS(nu,4) / LTAU(nu,5))
     +                 * dsin(LW(nu) * LTS(nu,4) + LTHETA(nu)))
     +                 / dexp(LTS(nu,4) / LTAU(nu,4))
            LSIG(nu,2) = dsqrt(abs((LTS(nu,5) - LTS(nu,6))**2/(-2.d0)
     +                 / dlog((LCO(nu,5)*dexp(LTS(nu,5) / LTAU(nu,4))
     +                 + LCO(nu,6) * dexp(LTS(nu,5)
     +                 / LTAU(nu,5))*dsin(LW(nu)*LTS(nu,5)
     +                 + LTHETA(nu))) / LCO(nu,7))))
            LSIG(nu,3) = LSIG(nu,2) * LSIG(nu,3)
            IF (LSIG(nu,2) .NE. 0) THEN
              LCO(nu,9) = LCO(nu,7)
     +                  * dexp((LTS(nu,7) - LTS(nu,6))**2
     +                  * (-0.5d0) / LSIG(nu,3)**2)
     +                  / dexp(LTS(nu,7) / LTAU(nu,6))
            ELSE
              LCO(nu,9) = (LCO(nu,5)
     +                  * dexp(LTS(nu,5) / LTAU(nu,4))
     +                  + LCO(nu,6)
     +                  * dexp(LTS(nu,5) / LTAU(nu,5))
     +                  * dsin(LW(nu) * LTS(nu,5)
     +                  + LTHETA(nu))) / dexp(LTS(nu,7) / LTAU(nu,6))
            END IF
            LCO(nu,8) = dsqrt(LCO(nu,7)**2
     +                - (LCO(nu,9)
     +                * dexp(LTS(nu,6) / LTAU(nu,6)))**2
     +                - (LCO(nu,10)
     +                / ((1.d0 + LTS(nu,6) / LDELT(nu))**LTN(nu)))**2)
            LCO(nu,11) = LCO(nu,10)
     +                 / ((1.d0 + LTS(nu,9) / LDELT(nu))**LTN(nu))

* determine needed constants for neutrino temperature function

            BSIGB(nu) = BSIGA(nu)*BSIGB(nu)
            BC(nu,5)  = BC(nu,4) + BC(nu,2) + BC(nu,3)*BTB(nu)
     +                - BC(nu,6)/(1.d0 + BTB(nu)/BDT(nu))**BN(nu)

         END DO  ! end loop over 3 groups of neutrinos
         CLOSE (unit = 13)

c         fermi_file = 'dat/model/Burrows/fermi_Burrows_eta.dat'

*****************************************************************************
* Beacom & Vogel
*****************************************************************************
      ELSEIF (sn_model_switch .EQ. 2) THEN
* energy (MeV)
         emin  = 1.00000d-01
         emax  = 1.70000d+02
         estep = 1.00000d-01
         energy_entries = nint((emax - emin) / estep + 1)
* time (ms)
c         tstep1      = 5.00000d+00
         tstep1      = 10.00000d+00  ! for Mac OS X
         time_start1 = -9.00000d+01
         time_stop1  = 2.00000d+04
         tstep2      = 5.00000d+00   ! need for sn_total_energy
         time_start2 = 0.
         time_stop2  = 2.00000d+04   ! need for sampling limits
         time_entries1 = 0           ! not used here
         time_entries2 = 0           ! not used here
         time_entries = nint((time_stop1 - time_start1) / tstep1 + 1)

* read BeacomVogel model parameters

         OPEN (unit = 13, 
     +         file = 'dat/model/BeacomVogel/BeacomVogel.dat',
     +         status = 'old')

         READ (13,*) description, E_bind
         READ (13,*) description, tau1
         READ (13,*) description, tau2
         READ (13,*) description, temp_nu_e
         READ (13,*) description, temp_nu_eb
         READ (13,*) description, temp_nu_mu

         CLOSE (unit = 13)

*****************************************************************************
* Bruenn_15
*****************************************************************************
      ELSEIF (sn_model_switch .EQ. 3) THEN
* energy (MeV)
         emin  = 3.5435d+00
         emax  = 3.5435d+02
         estep = 0.1d0        ! not meaningful; needed in nu energy sampling
         energy_entries = file_energy_size
* time (ms)
         time_start1 = 1.0216d+00
         time_stop1  = 5.00000d+02
         tstep1      = 1.00000d+00
         time_start2 = 6.00000d+02
         time_stop2  = 3.4910d+03
         tstep2      = 1.00000d+02
         time_entries1 = 0
         time_entries2 = 0
         time_entries = file_time_size

*****************************************************************************
* Bruenn_25
*****************************************************************************
      ELSEIF (sn_model_switch .EQ. 4) THEN
* energy (MeV)
         emin  = 3.5435d+00
         emax  = 3.5435d+02
         estep = 0.1d0        ! not meaningful; needed in nu energy sampling
         energy_entries = file_energy_size 
* time (ms)
         time_start1 = 1.0249d+00
         time_stop1  = 5.00000d+02
         tstep1      = 1.00000d+00
         time_start2 = 6.00000d+02
         time_stop2  = 9.9100d+02
         tstep2      = 1.00000d+02
         time_entries1 = 0
         time_entries2 = 0
         time_entries = file_time_size

*****************************************************************************
* Mezzacappa 13
*****************************************************************************
      ELSEIF (sn_model_switch .EQ. 5) THEN
* energy (MeV)
         emin  = 5.d+00
         emax  = 3.d+02
         estep = 0.1d0        ! not meaningful; needed in nu energy sampling
         energy_entries = file_energy_size
* time (ms)
         time_start1 = -8.995d+01
         time_stop1  = 5.00000d+02
         tstep1      = 1.00000d+00
         time_start2 = 6.00000d+02
         time_stop2  = 4.99609237d+02
         tstep2      = 1.00000d+02
         time_entries1 = 0
         time_entries2 = 0
         time_entries = file_time_size
      ENDIF  ! end models

C  Dummy energy variables  (integer values for DO loops)

      dum_emin = nint(emin / estep)
      dum_emax = nint(emax / estep)
      dum_estep = 1
     
C  Dummy time variables (integer values for DO loops)

      dum_tmin1 = nint(time_start1 / tstep1)
      dum_tmax1 = nint(time_stop1 / tstep1)
      dum_tstep1 = 1

      dum_tmin2 = nint(time_start2 / tstep2)
      dum_tmax2 = nint(time_stop2 / tstep2)
      dum_tstep2 = 1

C  Supernova Flux Entries

      IF (sn_model_switch .EQ. 1) THEN     ! Burrows
         flux_entries1 = time_entries1 * energy_entries
         flux_entries2 = time_entries2 * energy_entries
         flux_entries = flux_entries1 + flux_entries2
      ELSEIF (sn_model_switch .EQ. 2) THEN ! Beacom & Vogel
         flux_entries = time_entries * energy_entries
      ELSEIF (sn_model_switch .GT. 2) THEN
         flux_entries = file_flux_size
      ENDIF

C  Number of Supernova to Generate

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) n_supernova

C  Flavour Swapping Flag

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) swap_flag

C  Supernova Generator Flag

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) gen_flag
      READ (10,9000) description
      READ (10,1002) counts_file

C  Supernova Luminosity/Energy Fit Output Flag and Filename
      
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) lumen_fit_file_flag
      READ (10,9000) description
      READ (10,1002) lumen_nu_e_file
      READ (10,9000) description
      READ (10,1002) lumen_nu_eb_file
      READ (10,9000) description
      READ (10,1002) lumen_nu_mu_file

C  Supernova Luminosity Cutoff Flag and Time

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) lum_cutoff_flag
      READ (10,9000) description
      READ (10,1000) lum_cutoff_time

C  Supernova Fermi-Dirac Shape File
      
      READ (10,9000) description
      READ (10,9000) description  
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) FD_shape_file_flag
      READ (10,9000) description
      READ (10,1000) FD_time
      READ (10,9000) description   
      READ (10,1002) FD_file_nu_e
      READ (10,9000) description
      READ (10,1002) FD_file_nu_eb
      READ (10,9000) description   
      READ (10,1002) FD_file_nu_mu

C  Supernova Central Energy

      READ (10,9000) description
      READ (10,9000) description       
      READ (10,9000) description      
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description     
      READ (10,1001) central_flag
      READ (10,9000) description
      READ (10,1002) central_file

C  Supernova Cross Section File

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      
      READ (10,9000) description
      READ (10,1001) xsection_file_flag
      READ (10,9000) description
      READ (10,1002) xsection_file

C  Detector Elastic Scattering Electron Energy Threshold

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1000) T_e_min

C  Supernova Neutrino Countrate Sum Files
      
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description 
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1002) all_crfile   
      READ (10,9000) description
      READ (10,*) es_crfile, cc_crfile, nc_crfile
      READ (10,9000) description
      READ (10,*) nu_e_crfile, nu_eb_crfile, nu_mu_crfile   
      READ (10,9000) description
      READ (10,*) d2o_crfile, h2o_crfile, av_crfile
      
C  Supernova Neutrino Count Sum Files

      READ (10,9000) description
      READ (10,9000) description
      READ (10,1002) all_cfile
      READ (10,9000) description
      READ (10,*) es_cfile, cc_cfile, nc_cfile
      READ (10,9000) description
      READ (10,*) nu_e_cfile, nu_eb_cfile, nu_mu_cfile
      READ (10,9000) description
      READ (10,*) d2o_cfile, h2o_cfile, av_cfile

C  Supernova Neutrino Cumulative Count Sum Files

      READ (10,9000) description
      READ (10,9000) description
      READ (10,1002) all_cumcfile
      READ (10,9000) description
      READ (10,*) es_cumcfile, cc_cumcfile, nc_cumcfile
      READ (10,9000) description
      READ (10,*) nu_e_cumcfile, nu_eb_cumcfile, nu_mu_cumcfile
      READ (10,9000) description
      READ (10,*) d2o_cumcfile, h2o_cumcfile, av_cumcfile

C  Poisson Statistics Flag

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description 
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) pois_flag

C  DAQ Trigger Gate Stuff (trigger_gate in [ms])

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) time_overlap_flag
      READ (10,9000) description
      READ (10,1002) time_overlap_file
      READ (10,9000) description
      READ (10,1000) trigger_gate

C  Supernova Sampling Files

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description   
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1002) nu_sample_file
      READ (10,9000) description
      READ (10,1002) part_sample_file

C  Supernova Sampling Efficiency Flag

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) efficiency_flag

C  SNOMAN "MCPL" output file stuff

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) mcpl_flag
      READ (10,9000) description
      READ (10,1002) mcpl_file
      READ (10,9000) description
      READ (10,1002) mcpl_default

C  SNOMAN ".cmd" output file stuff

      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1001) cmd_flag
      READ (10,9000) description
      READ (10,1002) cmd_default
      READ (10,9000) description
      READ (10,1001) zdab_flag
      READ (10,9000) description
      READ (10,1002) zdab_default
      READ (10,9000) description
      READ (10,1002) hbook_default

C  SNO NHIT Threshold
      
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description
      READ (10,9000) description

      READ (10,9000) description
      READ (10,1000) nhit_thres
      READ (10,9000) description
      READ (10,1000) nhit_thres
      READ (10,9000) description
      READ (10,1000) pmt_efficiency
      READ (10,9000) description
      READ (10,1000) h2o_in_d2o
      READ (10,9000) description
      READ (10,1000) O17_in_d2o

**************************************************************************
*   Formats for READ statements
**************************************************************************

1000  FORMAT (1x,e16.9)
1001  FORMAT (1x,i10)
1002  FORMAT (1x,a40)
1003  FORMAT (1x,a80)
1099  FORMAT (1x,i4,i2,i2,1x,i2,i2,e16.9)
9000  FORMAT (a80)

      CLOSE (unit = 10)

      RETURN
      END
