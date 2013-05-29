***********************************************************************
*  Filename: sn_flux.f
***********************************************************************
*
*   Calculates the normalized flux of neutrinos at the SN source. 
*   The model parameters from various authors can be found in the 
*   "dat/model" area (ie., luminosity models, neutrino temperatures, etc.)
*
*    Inputs: 
*	o  this interaction id (not array)
*	o  neutrino energy  [MeV]
*	o  supernova (undelayed) time  [ms]
*
*    Outputs: 
*	o  supernova source flux  [(#/(s*MeV)]
*
***********************************************************************

      FUNCTION sn_flux(this_interaction,E_nu,time)

      IMPLICIT none

      INCLUDE 'flux_size.inc'
      INCLUDE 'file_flux_size.inc'
      INCLUDE 'id_nu.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON           /sn_model/ 
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /sn_flux_entries/
     +                 flux_entries1, flux_entries2, flux_entries
      INTEGER          flux_entries1, flux_entries2, flux_entries

      COMMON           /sn_FD_shape/ 
     +                 FD_shape_file_flag,
     +                 FD_nu_e_burst_file, FD_nu_e_cooling_file,
     +                 FD_nu_eb_burst_file, FD_nu_eb_cooling_file,
     +                 FD_nu_mu_burst_file, FD_nu_mu_cooling_file
      INTEGER          FD_shape_file_flag
      CHARACTER*40     FD_nu_e_burst_file, FD_nu_e_cooling_file,
     +                 FD_nu_eb_burst_file, FD_nu_eb_cooling_file,
     +                 FD_nu_mu_burst_file, FD_nu_mu_cooling_file

      COMMON           /L_file/
     +                 Le, Leb, Lmu
      REAL*8           Le(file_flux_size),
     +                 Leb(file_flux_size),
     +                 Lmu(file_flux_size)
      
      INTEGER          swap_flag
      COMMON           swap_flag

C  Global Variable Declarations

      INTEGER          this_interaction        ! input
      REAL*8           E_nu, time              ! input

C  Local Variable Declarations

c      REAL*8           time_delay, time_d ! don't use delayed time for now
       REAL*8           sn_flux, sn_num, fd_spectrum         ! functions
       REAL*8           sn_flux_nu(flux_size)       ! Burrows,Beacom&Vogel
c      REAL*8           sn_flux_nu(file_flux_size) ! Bruenn, Mezzacappa

      LOGICAL          same_interaction
      INTEGER          nu_type
      INTEGER          nu_type_swap,nu_type_old
      INTEGER          flux_index

      INTEGER          n_calls,
     +                 last_interaction, last_nu, last_int_nu

      SAVE             n_calls,
     +                 last_interaction, last_nu, last_int_nu
 
*************************************************************************
*  Initializations
*************************************************************************

      DATA             same_interaction /.false./,
     +                 n_calls          /0/,
     +                 last_interaction /0/,
     +                 last_nu          /0/ 

      DATA             sn_flux_nu /flux_size * 0.d0/  ! Burrows, Beacom&Vogel
c      DATA             sn_flux_nu /file_flux_size * 0.d0/ ! Bruenn, Mezzacappa

      IF (this_interaction .EQ. last_interaction) THEN
         same_interaction = .true.
      ELSE
         same_interaction = .false.
         n_calls = 0                ! zero counter for each interaction
         last_int_nu = last_nu      ! set last interaction neutrino type
      ENDIF

*************************************************************************

      n_calls = n_calls + 1

*************************************************************************
* If we have already calculated the sn_flux once, then other interactions
* which have the same nu_type will have the same sn_flux
*************************************************************************

* get neutrino type so we use the appropriate supernova flux array

      CALL get_nu_type(this_interaction,nu_type)

* determine position in flux array for given neutrino energy and time
* (only one neutrino type is considered at a time, so flux_index does
* not have to be the same for all neutrino types)

c      WRITE(*,*) E_nu, time

      CALL get_flux_index(E_nu,time,flux_index)

* if done integration, use stored flux array for sampling

      IF (      same_interaction
     +    .AND. n_calls .GT. flux_entries) THEN
         IF (n_calls .EQ. flux_entries + 1)
     +    WRITE (*,*) '  --> using stored array sn_flux value for sampli
     +ng'
         sn_flux = sn_flux_nu(flux_index)

* if new interaction with same neutrino type, use stored flux array

      ELSEIF (.NOT. same_interaction
     +        .AND. nu_type .EQ. last_int_nu) THEN
c         IF (n_calls .EQ. 1)
c     +    WRITE (*,*) '  --> using stored array sn_flux value'
         sn_flux = sn_flux_nu(flux_index)

* need to perform calculation

      ELSE
c         IF (n_calls .EQ. 1) THEN
c            IF (nu_type .EQ. id_nu_e) THEN
c               WRITE (*,*) '  --> calculating SN src flux for electron n
c     +eutrinos'
c            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
c               WRITE (*,*) '  --> calculating SN src flux for electron a
c     +ntineutrinos'
c            ELSEIF (nu_type .EQ. id_nu_mu) THEN
c               WRITE (*,*) '  --> calculating SN src flux for muon neutr
c     +inos'
c            ELSEIF (nu_type .EQ. id_nu_mubar) THEN
c               WRITE (*,*) '  --> calculating SN src flux for muon antin
c     +eutrinos'
c            ELSEIF (nu_type .EQ. id_nu_tau) THEN
c               WRITE (*,*) '  --> calculating SN src flux for tau neutri
c     +nos'
c            ELSEIF (nu_type .EQ. id_nu_taubar) THEN
c               WRITE (*,*) '  --> calculating SN src flux for tau antine
c     +utrinos'
c            ENDIF
c         ENDIF


*************************************************************************
* Preform flavour swapping if swap_flag = 1
*************************************************************************

         nu_type_old = nu_type
      IF (swap_flag .EQ. 1) THEN
         CALL flavour_swap(nu_type,E_nu,time,nu_type_swap)
         nu_type = nu_type_swap
      ENDIF

*************************************************************************
* Burrows et al. (1992) generic model or Beacom & Vogel generic model
*************************************************************************

         IF (     sn_model_switch .EQ. 1   ! Burrows
     +       .OR. sn_model_switch .EQ. 2)  ! Beacom & Vogel
     +   THEN

c* if massive neutrinos then adjust for energy-dependent time delay;
c* time_delay returns undelayed time if no massive neutrinos
c
c            time_d = time_delay(this_interaction,E_nu,time)

***********************************************************************
*  Calculate full supernova neutrino flux, combining the number of
*  neutrinos and the Fermi-Dirac spectrum
*    o sn_num --> numbers of neutrinos for the neutrino species
*    o fermi_spectrum --> Fermi-Dirac supernova neutrino spectrum for the
*                         neutrino species
*  Note: number and spectrum calculations are done using "star time" (ie., 
*        undelayed time)
***********************************************************************

C       WRITE (*,*) sn_num(nu_type,time)
C       WRITE (*,*) fd_spectrum(this_interaction,nu_type,E_nu,time)


            sn_flux =  sn_num(nu_type,time)
     +               * fd_spectrum(this_interaction,nu_type,E_nu,time)

* fill sn_flux array since we haven't gone through the entire
* energy and time spectrum for one entire interaction yet
* [Note: we are explicitly assuming that the integrations over time and
*        energy are being performed the first time this routine is called]

* flux_index needs time_d?????? (not sure if this distribution is sampled
* correctly; for now just keep undelayed references to array indices)

         IF (swap_flag .EQ. 1 .AND. 
     +       nu_type_old .NE. nu_type_swap) THEN
            IF (nu_type .EQ. id_nu_e .OR.
     +          nu_type .EQ. id_nu_ebar) THEN
                sn_flux = 0.5*sn_flux
            ELSE
                sn_flux = 2*sn_flux
            ENDIF
         ENDIF

            sn_flux_nu(flux_index) = sn_flux

*************************************************************************
**************************************************************************
*  Bruenn (2000) models (include Mezzacappa for now)
**************************************************************************

         ELSEIF (     sn_model_switch .EQ. 3
     +           .OR. sn_model_switch .EQ. 4
     +           .OR. sn_model_switch .EQ. 5) THEN

* assign fluxes based on the neutrino type of the current interaction

            IF (nu_type .EQ. id_nu_e) THEN
               sn_flux = Le(flux_index)

c               WRITE (*,*) 'sn_flux =', sn_flux

            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
               sn_flux = Leb(flux_index)

            ELSEIF (     nu_type .EQ. id_nu_mu
     +              .OR. nu_type .EQ. id_nu_mubar
     +              .OR. nu_type .EQ. id_nu_tau
     +              .OR. nu_type .EQ. id_nu_taubar) THEN
               sn_flux = Lmu(flux_index)
            ENDIF

* fill sn_flux array since we haven't gone through the entire
* energy and time spectrum for one entire interaction yet

         IF (swap_flag .EQ. 1 .AND. 
     +       nu_type_old .NE. nu_type_swap) THEN
            IF (nu_type .EQ. id_nu_e .OR. 
     +          nu_type .EQ. id_nu_ebar) THEN
                sn_flux = 0.5*sn_flux
            ELSE
                sn_flux = 2*sn_flux
            ENDIF
         ENDIF

           sn_flux_nu(flux_index) = sn_flux

         ENDIF  ! models
      ENDIF  ! test whether we can use the last interaction sn_flux

* remember some things about this time through the routine

      last_interaction = this_interaction
      last_nu = nu_type

c      WRITE(*,*) sn_flux, E_nu, time

      END
