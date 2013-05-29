***********************************************************************
*  Filename: fd_spectrum.f
***********************************************************************
*
*  This function subprogram assigns temperature and "chemical potential" 
*  parameters for the various neutrino types and calculates the Fermi-Dirac 
*  spectrum shape.  The normalization for the spectrum shape is calculated
*  by the "fd_norm" function.
*
*  Inputs:
*	o  this interaction id number (not array)
*	o  neutrino energy (MeV)
*       o  neutrino emission time (ms)
*       o  nu_type
*  Outputs:
*       o  Fermi-Dirac spectrum for the given neutrino energy and 
*          emission time
*       o  optionally, also writes the spectral values to a file specified
*          by the user
*
***********************************************************************
     
      FUNCTION fd_spectrum(this_interaction,nu_type,E_nu,this_time)

      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /Burrows_eta_stuff/
     +                 a, to, b
      REAL*8           a(3), to(3), b(3)

      COMMON           /Burrows_temp_stuff/
     +                 c, ta, tb, siga, sigb, n, dt
      REAL*8           c(3,6), ta(3),
     +                 tb(3), siga(3), sigb(3), n(3), dt(3)

      COMMON           /BeacomVogel_temp_stuff/
     +                 temp_nu_e, temp_nu_eb, temp_nu_mu
      REAL*8           temp_nu_e, temp_nu_eb, temp_nu_mu

      COMMON           /FD_shape/
     +                 FD_time, FD_shape_file_flag,
     +                 FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu
      INTEGER          FD_shape_file_flag
      REAL*8           FD_time
      CHARACTER*40     FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu

      COMMON           /sn_energy_entries/
     +                 energy_entries
      INTEGER          energy_entries

C  Global Variable Declarations

      INTEGER          this_interaction          ! input
      INTEGER          nu_type                   ! input
      REAL*8           E_nu, this_time           ! input

C  Local Variable Declarations

      REAL*8           fd_spectrum               ! function

      INTEGER          nu                        ! Burrows
      INTEGER          last_nu_type
      REAL*8           t                         ! Burrows
c      INTEGER          j                        ! (Burrows integration)
c      REAL*8           sum                      ! (Burrows integration)

      REAL*8           temp_nu, eta_nu, norm, fd_norm

      INTEGER          E_entries

      INTEGER          calls_nu_e, calls_nu_eb, calls_nu_mu,
     +                 last_interaction
      REAL*8           last_time, last_temp, last_eta, last_norm
      LOGICAL          first_call, done_nu_e, done_nu_eb, done_nu_mu

      SAVE             first_call, 
     +                 calls_nu_e, calls_nu_eb, calls_nu_mu, E_entries,
     +                 last_interaction, last_time, last_nu_type,
     +                 last_temp, last_eta, last_norm,
     +                 done_nu_e, done_nu_eb, done_nu_mu

************************************************************************
*  Initializations
************************************************************************

      DATA             calls_nu_e        /0/,
     +                 calls_nu_eb       /0/,
     +                 calls_nu_mu       /0/,
     +                 E_entries         /0/,
     +                 last_interaction  /0/,
     +                 last_time         /0.d0/,
     +                 last_temp         /0.d0/,
     +                 last_eta          /0.d0/,
     +                 last_norm         /0.d0/

      DATA             first_call        /.true./,
     +                 done_nu_e         /.false./,
     +                 done_nu_eb        /.false./,
     +                 done_nu_mu        /.false./

************************************************************************

      IF (      this_interaction .EQ. last_interaction
     +    .AND. this_time .EQ. last_time .AND. nu_type .EQ.
     +    last_nu_type) THEN
         temp_nu = last_temp
         eta_nu  = last_eta
         norm    = last_norm
      ELSE

************************************************************************
* Burrows
************************************************************************

         IF (sn_model_switch .EQ. 1) THEN

* Burrows starts infall at time=0 rather than -100 ms and uses seconds;
* we'll convert to that convention
            t = (this_time + 100.d0) / 1000.d0

            nu = 0

            IF (nu_type .EQ. id_nu_e) THEN
               nu = 1
            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
               nu = 2
            ELSEIF (     nu_type .EQ. id_nu_mu
     +              .or. nu_type .EQ. id_nu_mubar
     +              .or. nu_type .EQ. id_nu_tau
     +              .or. nu_type .EQ. id_nu_taubar) THEN
               nu = 3
            ENDIF

           IF (t .LT. 0.d0) t = 0.d0

* assign eta values (depends only on time)

           IF (t .LT. to(nu)) THEN
             eta_nu = a(nu)
           ELSE
             eta_nu = b(nu)
           END IF

* assign neutrino temperatures (depends only on time and neutrino type)

           IF (t .LT. ta(nu)) THEN
             temp_nu = c(nu,1)
           ELSEIF (t .LT. tb(nu)) THEN
             temp_nu = c(nu,2) + c(nu,3)*t
     +             + c(nu,4)*dexp(((t - tb(nu))/siga(nu))**2/(-2.d0))
           ELSE
             temp_nu = c(nu,5)*dexp(((t - tb(nu))/sigb(nu))**2/(-2.d0))
     +                + c(nu,6)/(1.d0 + t/dt(nu))**n(nu)
           END IF

c           norm = 0.25d0*(eta_nu**4) + 4.9348022d0*(eta_nu**2) 
c     +           + 11.36439d0
c           SUM  = 0.d0
c           DO j = 1,5
c             sum = sum + (-1.d0)**j/dble(j**4) * dexp((-1.d0)*j*eta_nu)
c           END DO
c           norm = (norm + 6.d0*sum) * (temp_nu)**4

c           write (*,*) 'time, T, eta = ', this_time, temp_nu, eta_nu

*************************************************************************
* Beacom & Vogel model
*************************************************************************

         ELSEIF (sn_model_switch .EQ. 2) THEN

            IF (nu_type .EQ. id_nu_e) THEN
               temp_nu = temp_nu_e
            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
               temp_nu = temp_nu_eb
            ELSEIF (     nu_type .EQ. id_nu_mu
     +              .or. nu_type .EQ. id_nu_mubar
     +              .or. nu_type .EQ. id_nu_tau
     +              .or. nu_type .EQ. id_nu_taubar) THEN
               temp_nu = temp_nu_mu
            ENDIF

            eta_nu = 0.d0  ! no chemical potential in this model

         ENDIF   ! check for which SN model

* calculate Fermi-Dirac spectrum normalization

         norm = fd_norm(temp_nu,eta_nu)

* remember some things about this time through the routine

         last_interaction = this_interaction
         last_time = this_time
         last_temp = temp_nu
         last_eta  = eta_nu
         last_norm = norm
         last_nu_type = nu_type

      ENDIF  ! check if same interaction and same time (for normalization)

      fd_spectrum = (1.d0 / norm)
     +             * E_nu**2 / (1.d0 + dexp(E_nu / temp_nu - eta_nu))

**************************************************************************
*  Write the various Fermi-Dirac spectra to user specified file.  
*  -> need to check all energies and all times
*  (assumes time integration in "count_sum" done first, then energy)
**************************************************************************

      IF (FD_shape_file_flag .EQ. 1) THEN

* don't need to write these values again after we've written all of the
* energy entries already
*  -> Burrows: both the burst and the cooling (hence "2 * energy_entries")
*  -> Beacom & Vogel: just energy_entries

         IF (first_call) THEN
            first_call = .false.
            IF (sn_model_switch .EQ. 1) THEN      ! Burrows
               E_entries = 2 * energy_entries         
            ELSEIF (sn_model_switch .EQ. 2) THEN  ! Beacom & Vogel
               E_entries = energy_entries             
            ENDIF  ! check for model
         ENDIF     ! check for first call

C         CALL get_nu_type(this_interaction,nu_type)
  
* nu_e
         IF (.NOT. done_nu_e
     +       .AND. nu_type .EQ. id_nu_e) THEN

            IF (this_time .EQ. FD_time) THEN
               calls_nu_e = calls_nu_e + 1
               IF (calls_nu_e .EQ. 1) 
     +         WRITE (*,*) 'nu_e:temp_nu, eta_nu =', temp_nu, eta_nu
               WRITE (70,*) E_nu, fd_spectrum
               IF (calls_nu_e .EQ. E_entries) done_nu_e = .true.
            ENDIF

* nu_e_bar
         ELSEIF (.NOT. done_nu_eb
     +           .AND. nu_type .EQ. id_nu_ebar) THEN

            IF (this_time .EQ. FD_time) THEN
               calls_nu_eb = calls_nu_eb + 1
               IF (calls_nu_eb .EQ. 1) 
     +         WRITE (*,*) 'nu_eb:temp_nu, eta_nu =', temp_nu, eta_nu
               WRITE (71,*) E_nu, fd_spectrum
               IF (calls_nu_eb .EQ. E_entries) done_nu_eb = .true.
            ENDIF

*nu_mu
         ELSEIF (.NOT. done_nu_mu
     +           .AND. nu_type .EQ. id_nu_mu) THEN

            IF (this_time .EQ. FD_time) THEN
               calls_nu_mu = calls_nu_mu + 1
               IF (calls_nu_mu .EQ. 1) 
     +         WRITE (*,*) 'nu_mu:temp_nu, eta_nu =', temp_nu, eta_nu
               WRITE (72,*) E_nu, fd_spectrum
               IF (calls_nu_mu .EQ. E_entries) done_nu_mu = .true.
            ENDIF ! sort by time
         ENDIF    ! sort by neutrino type
      ENDIF       ! check to see if we write the spectra to file

      END
