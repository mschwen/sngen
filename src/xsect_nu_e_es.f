***********************************************************************
*  Filename: xsect_nu_e_es.f
***********************************************************************
*
*  This routine returns the total nu + e ES cross section for a given
*  neutrino energy.  Two sets of cross sections are included:
*	o Bahcall (1989) [complete with T_min dependence]
*	o Burrows et al. (1992)
*  Use the "author_params" flag to chose between the two.
*  Default is Bahcall.
*
*  Four cross sections types are included in this routine:
*	o nu_e, nu_e_bar, nu_"mu", nu_"mu_bar"
*  except in the case of Burrows et al. where they don't distinguish 
*  between nu_"mu" and nu_"mu_bar".
*
*  Note:  Bahcall's formulae are expressed in terms of the electron mass.
*     
*  Inputs:
*       o  interaction id  (not array)
*	o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy (m^2)
*
***********************************************************************
     
      SUBROUTINE xsect_nu_e_es(interaction,E_nu,xsect_total)

      IMPLICIT none

      INCLUDE 'energy_size.inc'
      INCLUDE 'xsect_size.inc'          ! new method
      INCLUDE 'aspline_size.inc'        ! new method
      INCLUDE 'id_nu.inc'
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON           /generator_flag/
     +                 gen_flag, counts_file
      INTEGER          gen_flag
      CHARACTER*40     counts_file

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /sn_model_params/
     +                 use_author_params
      INTEGER          use_author_params

      COMMON           /sn_xsection_file/
     +                 xsection_file_flag,xsection_file
      INTEGER          xsection_file_flag
      CHARACTER*40     xsection_file

      COMMON           /energy_limits/
     +                 emin, emax, estep,
     +                 dum_emin, dum_emax, dum_estep
      REAL*8           emin, emax, estep
      INTEGER          dum_emin, dum_emax, dum_estep

      COMMON           /sn_energy_entries/
     +                 energy_entries
      INTEGER          energy_entries

      COMMON           /energy_index/   ! new
     +                 this_energy_index
      INTEGER          this_energy_index
 
      COMMON           /det_threshold/
     +                 E_min
      REAL*8           E_min

C  Global Variable Declarations

      INTEGER          interaction  ! input
      REAL*8           E_nu         ! input
      REAL*8           xsect_total  ! output   

C  Local Variable Declarations

      REAL*8           energy, xsect(energy_size) 
      INTEGER          e_dum, e_index

      INTEGER          nu_type

* new method
      LOGICAL          first_call
      INTEGER          i

      INTEGER          xs_data_entries
      PARAMETER        (xs_data_entries = xs_nu_e_entries)        ! default
      CHARACTER*40     xs_data_file
      REAL*8           E_data(xs_data_entries), 
     +                 xs_nue_data(xs_data_entries), 
     +                 xs_numu_data(xs_data_entries)
      REAL*8           b_nue_xs(aspline_size), 
     +                 c_nue_xs(aspline_size),
     +                 d_nue_xs(aspline_size)
      REAL*8           b_numu_xs(aspline_size), 
     +                 c_numu_xs(aspline_size),
     +                 d_numu_xs(aspline_size)
      REAL*8           xsect_fit

      REAL*8           Eth                      ! Burrows
      DOUBLE PRECISION T_e_min, T_e_max
      DOUBLE PRECISION g_L, g_R, g_L2, g_R2
      DOUBLE PRECISION term1, term2, term3
      DOUBLE PRECISION sig_e

      INTEGER          n_calls, last_interaction
      
      SAVE             n_calls, last_interaction
      SAVE             xsect

********************************************************************
*  Initializations
********************************************************************
      
      DATA n_calls           /0/
      DATA last_interaction  /0/

      DATA E_data              /xs_data_entries * 0./,
     +     xs_nue_data         /xs_data_entries * 0./,
     +     xs_numu_data        /xs_data_entries * 0./
      
      DATA b_nue_xs            /aspline_size * 0./,
     +     c_nue_xs            /aspline_size * 0./,
     +     d_nue_xs            /aspline_size * 0./,
     +     b_numu_xs           /aspline_size * 0./,
     +     c_numu_xs           /aspline_size * 0./,
     +     d_numu_xs           /aspline_size * 0./

      DATA xsect               /energy_size * 0./

      xs_data_file = 'dat/xsect/xs_nu_e_es_rad.dat'           ! default

*********************************************************************

* reset counter for the next interaction

      IF (interaction .NE. last_interaction) n_calls = 0

* increment counter for this interaction

      n_calls = n_calls + 1

* if we've gone through the energy loop once, just use saved values
* (doesn't work if we bypass the integrations over time and energy by
*  reading in the count sum totals, therefore always calculate the cross
*  section explicitly in that case)
      
      IF (      (n_calls .GT. energy_entries)
     +    .AND. (gen_flag .NE. 3)            ) THEN

c         write (*,*) '  --> using saved xsect values!'
c         write (*,*) '      n_calls, energy_entries=',
c     +                n_calls, energy_entries

         xsect_total = xsect(this_energy_index)

c         write (*,*) 'this_energy_index, xsect_total=',
c     +                this_energy_index, xsect_total

* if we have filled the 'xsect' array, then write out the cross sections
* to file (once)

         IF (      xsection_file_flag .EQ. 1
     +       .AND. n_calls .EQ. energy_entries + 1) THEN

            CALL file_open(51,xsection_file,interaction)

            DO 10 e_dum = dum_emin, dum_emax, dum_estep
               energy = float(e_dum) * estep
               e_index = nint((energy - emin) / estep) + 1

               WRITE (51,*) energy, xsect(e_index)
10          CONTINUE

            CLOSE (unit = 51)
         ENDIF
         
* we haven't saved the values yet, so do the calculation

      ELSE

         CALL get_nu_type(interaction,nu_type)

*************************************************************************
*  Burrows cross section (courtesy Adam Burrows)
*  (g_L and g_R variable names are used for convenience)
*  Note: Burrows imposes a 5 MeV threshold on all SNO reactions
*************************************************************************

         IF (      sn_model_switch .EQ. 1
     +       .AND. use_author_params .EQ. 1) THEN

            T_e_min = 5.d0   ! 5 MeV detector threshold

            IF (nu_type .EQ. id_nu_e) THEN

               g_L = 2.1316d0
               g_R = 0.2116d0

            ELSEIF (nu_type .EQ. id_nu_ebar) THEN

               g_L = 0.2116d0
               g_R = 2.1316d0

* Burrows treats all the nu_"mu" elastic scattering cross sections the same
            ELSEIF (     nu_type .EQ. id_nu_mu
     +              .OR. nu_type .EQ. id_nu_mubar
     +              .OR. nu_type .EQ. id_nu_tau
     +              .OR. nu_type .EQ. id_nu_taubar) THEN

               g_L = 0.2516d0
               g_R = 0.2516d0 

            ENDIF

            xsect_total = 0.d0

            Eth = T_e_min / (E_nu + 0.000001d0)

            IF (Eth .LE. 1.d0)
     +       xsect_total = E_nu / 8.d0 * ((g_L + g_R / 3.d0)
     +                    - Eth * ((g_L + g_R)
     +                    - g_R * Eth * (1.d0 - Eth / 3.d0)))   

            xsect_total = 3.337d0 * xsect_total * 1.d-44 * cm2m**2

         ELSE  ! use the default method

            IF (      sn_model_switch .EQ. 2         ! Beacom & Vogel
     +          .AND. use_author_params .EQ. 1) THEN
               T_e_min = 5.d0                        ! detector threshold
            ELSE
               T_e_min = E_min
            ENDIF

************************************************************************
*  Bahcall's method:
*  Assign the coupling constants for the various neutrino types
************************************************************************

* radiative corrections
            IF (     nu_type .EQ. id_nu_e
     +          .OR. nu_type .EQ. id_nu_mu
     +          .OR. nu_type .EQ. id_nu_tau) THEN

* nu_e
c               g_L = (+0.5d0 + sin2theta)
c               g_R = sin2theta

* nu_mu, nu_tau
c               g_L = (-0.5d0 + sin2theta)
c               g_R = sin2theta

               IF (first_call) THEN

                  first_call = .false.
 
                  OPEN (unit = 50, file = xs_data_file, status = 'old')
      
                  DO 20 i = 1, xs_data_entries
                     READ (50,*) E_data(i), 
     +                           xs_nue_data(i), xs_numu_data(i)
            
                     xs_nue_data(i)    = dlog10(xs_nue_data(i))
                     xs_numu_data(i) = dlog10(xs_numu_data(i))

20                CONTINUE

                  CLOSE (unit = 50)

* nue
                  CALL aspline(xs_data_entries,E_data,xs_nue_data,
     +                         b_nue_xs,c_nue_xs,d_nue_xs)
               
* numu, nutau
                  CALL aspline(xs_data_entries,E_data,xs_numu_data,
     +                         b_numu_xs,c_numu_xs,d_numu_xs)

               ENDIF     ! condition for first call

*****************************************************************************
*  Evaluate the Akima spline fit (done for every neutrino for each
*  interaction)
*****************************************************************************
* nue 
               IF (     interaction .EQ. 1
     +             .OR. interaction .EQ. 2
     +             .OR. interaction .EQ. 57
     +             .OR. interaction .EQ. 78) THEN

                  CALL aspline_eval(xs_data_entries,
     +                              E_data,xs_nue_data,
     +                              b_nue_xs,c_nue_xs,d_nue_xs,
     +                              E_nu,xsect_fit)

* numu, nutau
               ELSEIF (     interaction .EQ. 21
     +                 .OR. interaction .EQ. 22
     +                 .OR. interaction .EQ. 39
     +                 .OR. interaction .EQ. 40
     +                 .OR. interaction .EQ. 66
     +                 .OR. interaction .EQ. 72
     +                 .OR. interaction .EQ. 88
     +                 .OR. interaction .EQ. 94) THEN

                  CALL aspline_eval(xs_data_entries,
     +                              E_data,xs_numu_data,
     +                              b_numu_xs,c_numu_xs,d_numu_xs,
     +                              E_nu,xsect_fit)

               ENDIF  ! check for which interaction

               xsect_total = 10.d0**(xsect_fit)

            ELSE      ! reactions that do not have radiative corrections

               IF (nu_type .EQ. id_nu_ebar) THEN

                  g_L = sin2theta
                  g_R = (+0.5d0 + sin2theta)

               ELSEIF (     nu_type .EQ. id_nu_mubar
     +                 .OR. nu_type .EQ. id_nu_taubar) THEN

                  g_L = sin2theta
                  g_R = (-0.5d0 + sin2theta)

               ENDIF

               g_L2 = g_L**2
               g_R2 = g_R**2

**************************************************************************
*  Cross section energy dependent stuff (including T_min dependence)
**************************************************************************

* maximum electron kinetic energy

               T_e_max = 2.d0 * E_nu**2 / (2.d0 * E_nu + m_electron)

* various cross section terms (see Bahcall)

               term1 = (g_L2 + g_R2) * (T_e_max - T_e_min)
               term2 = (g_R2 + (g_L * g_R * m_electron) 
     +                / (2.d0 * E_nu)) * 
     +                 (T_e_max**2 - T_e_min**2) / E_nu
               term3 = g_R2 / (3.d0 * E_nu**2) * 
     +                 (T_e_max**3 - T_e_min**3)

* overall cross section constant
* (numerical value = 1.72329852464e-44 for PDG 2000 constants)

               sig_e = 2.d0 * (G_fermi**2 * hbar_c2) * m_electron / pi

* add all the terms and factor in the overall cross section constant

               xsect_total = sig_e * (term1 - term2 + term3)

            ENDIF     ! check for nu_e or nu_mu (radiative corrections)

* convert cross section from [cm] to [m]

            xsect_total = xsect_total * cm2m**2

* using detector threshold can give negative values, so set to 0 in that case

            IF (xsect_total .LT. 0.d0) xsect_total = 0.d0

        ENDIF  ! check for use of author parameters

* fill cross section array since we haven't gone through the entire 
* energy spectrum yet

         xsect(this_energy_index) = xsect_total

C         WRITE (*,*) 'this_energy_index, E_nu, xsect_total =',
C     +                this_energy_index, E_nu, xsect_total

      ENDIF  ! use saved values or do calculation

* remember some things about this time through the routine

      last_interaction = interaction

      RETURN
      END
