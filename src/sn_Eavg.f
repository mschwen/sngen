***********************************************************************
*  Filename: sn_Eavg.f
***********************************************************************
*
*  Inputs:
*	o  neutrino type 
*	o  neutrino emission time (ms)
*  Outputs:
*       o SN average neutrino energy
*
*  Note: this routine only depends on time, so when integrating, we should
*  save the previous number 
* 
***********************************************************************
     
      FUNCTION sn_Eavg(nu_type,this_time)

      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'lumen_size.inc'
      INCLUDE 'aspline_size.inc'

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

C  Global Variable Declarations

      INTEGER          nu_type                ! input
      REAL*8           this_time              ! input

C  Local Variable Declarations

      REAL*8           sn_Eavg  ! function
      INTEGER          j

      INTEGER          nu                          ! Burrows
      REAL*8           t,                
     +                 eta, beta, fth, sum, ftw 
      REAL*8           nu_temp                     ! Beacom & Vogel only

      LOGICAL          first_call

      INTEGER          Eavg_data_entries
      CHARACTER*40     Eavg_data_file

      INTEGER          i

      REAL*8           t_data(data_size), 
     +                 Eavg_e_data(data_size), Eavg_eb_data(data_size),
     +                 Eavg_mu_data(data_size), Eavg_mub_data(data_size)
      REAL*8           b_Eavg_e(aspline_size), 
     +                 c_Eavg_e(aspline_size),
     +                 d_Eavg_e(aspline_size)
      REAL*8           b_Eavg_eb(aspline_size),          
     +                 c_Eavg_eb(aspline_size),
     +                 d_Eavg_eb(aspline_size)
      REAL*8           b_Eavg_mu(aspline_size),          
     +                 c_Eavg_mu(aspline_size),
     +                 d_Eavg_mu(aspline_size)
      REAL*8           b_Eavg_mub(aspline_size),
     +                 c_Eavg_mub(aspline_size),
     +                 d_Eavg_mub(aspline_size)

      REAL*8           Eavg_fit

      REAL*8           last_time, last_sn_Eavg
      INTEGER          last_nu

      SAVE             first_call
      SAVE             Eavg_data_entries,
     +                 t_data, 
     +                 Eavg_e_data, Eavg_eb_data, 
     +                 Eavg_mu_data, Eavg_mub_data,
     +                 b_Eavg_e, c_Eavg_e, d_Eavg_e,
     +                 b_Eavg_eb, c_Eavg_eb, d_Eavg_eb,
     +                 b_Eavg_mu, c_Eavg_mu, d_Eavg_mu,
     +                 b_Eavg_mub, c_Eavg_mub, d_Eavg_mub

      SAVE             last_nu, last_time, last_sn_Eavg

************************************************************************
*  Initializations
************************************************************************

      DATA             first_call       /.true./

      DATA             t_data           /data_size * 0./,
     +                 Eavg_e_data      /data_size * 0./,
     +                 Eavg_eb_data     /data_size * 0./,
     +                 Eavg_mu_data     /data_size * 0./,
     +                 Eavg_mub_data    /data_size * 0./
      
      DATA             b_Eavg_e         /aspline_size * 0./,
     +                 c_Eavg_e         /aspline_size * 0./,
     +                 d_Eavg_e         /aspline_size * 0./,
     +                 b_Eavg_eb        /aspline_size * 0./,
     +                 c_Eavg_eb        /aspline_size * 0./,
     +                 d_Eavg_eb        /aspline_size * 0./,
     +                 b_Eavg_mu        /aspline_size * 0./,
     +                 c_Eavg_mu        /aspline_size * 0./,
     +                 d_Eavg_mu        /aspline_size * 0./,
     +                 b_Eavg_mub       /aspline_size * 0./,
     +                 c_Eavg_mub       /aspline_size * 0./,
     +                 d_Eavg_mub       /aspline_size * 0./

      DATA             last_time         /-99999./,
     +                 last_nu           /0/,
     +                 last_sn_Eavg      /0./

**************************************************************************

      IF (      this_time .EQ. last_time
     +    .AND. nu_type   .EQ. last_nu  ) THEN
         sn_Eavg = last_sn_Eavg
         RETURN
      ENDIF

****************************************************************************
****************************************************************************
* Burrows et al. (1992)
****************************************************************************

      IF (sn_model_switch .EQ. 1) THEN

         nu = 0
         sn_Eavg = 0.d0

* Burrows starts infall at time=0 and uses seconds;
* we'll convert to that convention
         t   = (this_time + 100.d0) / 1000.d0

         IF (nu_type .EQ. id_nu_e) THEN
            nu = 1
         ELSEIF (nu_type .EQ. id_nu_ebar) THEN
            nu = 2
         ELSEIF (     nu_type .EQ. id_nu_mu
     +           .or. nu_type .EQ. id_nu_mubar
     +           .or. nu_type .EQ. id_nu_tau
     +           .or. nu_type .EQ. id_nu_taubar) THEN
            nu = 3
         ENDIF

* assign eta values depending on time

        IF (T .LT. TO(NU)) THEN
          eta = A(NU)
        ELSE
          eta = B(NU)
        END IF

* beta

        IF (T .LT. TA(NU)) THEN
          BETA = C(NU,1)
        ELSE IF (T .LT. TB(NU)) THEN
          BETA = C(NU,2) + C(NU,3)*T
     $         + C(NU,4)*dexp(((T - TB(NU))/SIGA(NU))**2/(-2.d0))
        ELSE
          BETA = C(NU,5)*dexp(((T - TB(NU))/SIGB(NU))**2/(-2.d0))
     $         + C(NU,6)/(1.d0 + T/DT(NU))**N(NU)
        END IF
        BETA = 1.0/BETA

* functional form of the average energy curves (courtesey Adam Burrows)

        fth = 0.25d0*(eta**4) + 4.9348022d0*(eta**2) + 11.36439d0
        sum  = 0
        DO j = 1, 10
          sum = sum + (-1.d0)**j / dble(j**4) * dexp((-1.d0)*j*eta)
        END DO
        fth = fth + 6.d0 * sum
        ftw = eta**3/3.d0 + 3.290d0*eta
        sum = 0
        DO j = 1,10
          sum = sum + (-1.d0)**j / dble(j**3) * dexp((-1.d0)*j*eta)
        END DO
        ftw = ftw - 2.d0 * sum

        sn_Eavg = fth / ftw / beta

c        write (*,*) nu_type, this_time, nu_temp, sn_Eavg / nu_temp

*************************************************************************
* Beacom & Vogel model
*************************************************************************

      ELSEIF (sn_model_switch .EQ. 2) THEN

* average neutrino energy (depends on average neutrino temperature)

         IF (nu_type .EQ. id_nu_e) THEN
            nu_temp = temp_nu_e
         ELSEIF (nu_type .EQ. id_nu_ebar) THEN
            nu_temp = temp_nu_eb
         ELSEIF (     nu_type .EQ. id_nu_mu
     +           .OR. nu_type .EQ. id_nu_mubar
     +           .OR. nu_type .EQ. id_nu_tau
     +           .OR. nu_type .EQ. id_nu_taubar) THEN
            nu_temp = temp_nu_mu
         ENDIF

         sn_Eavg = 3.15d0 * nu_temp  ! 3.15 comes from FD statistics, no eta

*************************************************************************
* Detailed models
*************************************************************************

      ELSEIF (sn_model_switch .GT. 2) THEN

         IF (first_call) THEN

            first_call = .false.
 
            IF (sn_model_switch .EQ. 3) THEN

               Eavg_data_file = 
     +          'dat/model/Bruenn/15Ms/energy.dat'
               Eavg_data_entries = Eavg_Bruenn15_entries

            ELSEIF (sn_model_switch .EQ. 4) THEN

               Eavg_data_file = 
     +          'dat/model/Bruenn/25Ms/energy.dat'
               Eavg_data_entries = Eavg_Bruenn25_entries

            ELSEIF (sn_model_switch .EQ. 5) THEN

               Eavg_data_file = 
     +          'dat/model/Mezzacappa/13Ms/energy.dat'
               Eavg_data_entries = Eavg_Mezz13_entries

            ENDIF

            OPEN (unit = 50, file = Eavg_data_file, status = 'old')
      
            DO 20 i = 1, Eavg_data_entries
               READ (50,*) t_data(i), 
     +                     Eavg_e_data(i), Eavg_eb_data(i), 
     +                     Eavg_mu_data(i), Eavg_mub_data(i)

c               IF (t_data(i) - t_data(i-1) .LE. 0) THEN
c                  write (*,*) 'unordered!: i, t_data(i) =',i, t_data(i)
c               ENDIF

20          CONTINUE

            CLOSE (unit = 50)

* nu_e luminosity
            CALL aspline(Eavg_data_entries,t_data,Eavg_e_data,
     +                   b_Eavg_e,c_Eavg_e,d_Eavg_e)
 
* nu_eb luminosity
            CALL aspline(Eavg_data_entries,t_data,Eavg_eb_data,
     +                   b_Eavg_eb,c_Eavg_eb,d_Eavg_eb)

* nu_mu luminosity
            CALL aspline(Eavg_data_entries,t_data,Eavg_mu_data,
     +                   b_Eavg_mu,c_Eavg_mu,d_Eavg_mu)

* nu_mub luminosity
            CALL aspline(Eavg_data_entries,t_data,Eavg_mub_data,
     +                   b_Eavg_mub,c_Eavg_mub,d_Eavg_mub)

         ENDIF     ! condition for first call

* evaluate the splines for each call to the routine

         IF (nu_type .EQ. id_nu_e) THEN

            CALL aspline_eval(Eavg_data_entries,
     +                        t_data,Eavg_e_data,
     +                        b_Eavg_e,c_Eavg_e,d_Eavg_e,
     +                        this_time,Eavg_fit)

         ELSEIF (nu_type .EQ. id_nu_ebar) THEN

            CALL aspline_eval(Eavg_data_entries,
     +                        t_data,Eavg_eb_data,
     +                        b_Eavg_eb,c_Eavg_eb,d_Eavg_eb,
     +                        this_time,Eavg_fit)

         ELSEIF (     nu_type .EQ. id_nu_mu
     +           .OR. nu_type .EQ. id_nu_tau) THEN

            CALL aspline_eval(Eavg_data_entries,
     +                        t_data,Eavg_mu_data,
     +                        b_Eavg_mu,c_Eavg_mu,d_Eavg_mu,
     +                        this_time,Eavg_fit)

         ELSEIF (     nu_type .EQ. id_nu_mubar
     +           .OR. nu_type .EQ. id_nu_taubar) THEN

            CALL aspline_eval(Eavg_data_entries,
     +                        t_data,Eavg_mub_data,
     +                        b_Eavg_mub,c_Eavg_mub,d_Eavg_mub,
     +                        this_time,Eavg_fit)

         ENDIF   ! check for neutrino type

         sn_Eavg = Eavg_fit ! already in ergs/s

         IF (sn_Eavg .LE. 0.d0) THEN        ! some early Mezzacappa 
                                            ! fits are < 0; can't set = 0
            IF (nu_type .EQ. id_nu_e) THEN
               sn_Eavg = 5.6d0  ! modest changes have no effect
            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
               sn_Eavg = 6.1d0  ! modest changes have no effect
            ELSE
               sn_Eavg = 5.9d0  ! modest changes have no effect
            ENDIF
          ENDIF

      ENDIF   ! check for which model

* remember a few things about the last time this routine was called

      last_nu = nu_type
      last_time = this_time
      last_sn_Eavg = sn_Eavg

      END
