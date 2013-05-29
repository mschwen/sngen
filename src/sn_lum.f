***********************************************************************
*  Filename: sn_lum.f
***********************************************************************
*
* In all cases, "sn_lum" returns the luminosity for each of 
* 6 neutrinos types:
*   --> nu_e, nu_ebar, nu_mu, nu_mubar, nu_tau, nu_taubar
*
* That is, the "nu_mu" components are *not* combined into a single 
* luminosity value.
*
*  Inputs:
*	o  neutrino type 
*	o  neutrino emission time (ms)
*  Outputs:
*       o SN luminosity 
*
*  Note: this routine only depends on time, so when integrating, we should
*  save the previous number 
* 
***********************************************************************
     
      FUNCTION sn_lum(nu_type,this_time)

      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'convert_units.inc'
      INCLUDE 'lumen_size.inc'
      INCLUDE 'aspline_size.inc'

C  Common Block Declarations

      COMMON           /sn_model/ 
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /Burrows_lum_stuff/
     +                 ts, tau, co, sig, w, theta, tn, delt
      REAL*8           ts(3,0:9), tau(3,6), co(3,11), sig(3,4),
     +                 w(3), theta(3), tn(3), delt(3)

      COMMON           /BeacomVogel_lum_stuff/
     +                 E_bind, tau1, tau2
      REAL*8           E_bind, tau1, tau2

      COMMON           /lum_cutoff_stuff/
     +                 lum_cutoff_time,
     +                 lum_cutoff_flag
      INTEGER          lum_cutoff_flag
      REAL*8           lum_cutoff_time

C  Global Variable Declarations

      INTEGER          nu_type                ! input
      REAL*8           this_time              ! input

C  Local Variable Declarations

      REAL*8           sn_lum

      INTEGER          nu                ! Burrows neutrino type
      REAL*8           t                 ! Burrows time in seconds

      LOGICAL          first_call

      INTEGER          lum_data_entries
      CHARACTER*40     lum_data_file

      INTEGER          i

      REAL*8           t_data(data_size), 
     +                 lum_e_data(data_size), lum_eb_data(data_size),
     +                 lum_mu_data(data_size), lum_mub_data(data_size)
      REAL*8           b_lum_e(aspline_size), 
     +                 c_lum_e(aspline_size),
     +                 d_lum_e(aspline_size)
      REAL*8           b_lum_eb(aspline_size),          
     +                 c_lum_eb(aspline_size),
     +                 d_lum_eb(aspline_size)
      REAL*8           b_lum_mu(aspline_size),          
     +                 c_lum_mu(aspline_size),
     +                 d_lum_mu(aspline_size)
      REAL*8           b_lum_mub(aspline_size),
     +                 c_lum_mub(aspline_size),
     +                 d_lum_mub(aspline_size)

      REAL*8           lum_fit

      REAL*8           last_time, last_sn_lum
      INTEGER          last_nu

      SAVE             first_call
      SAVE             lum_data_entries,
     +                 t_data, 
     +                 lum_e_data, lum_eb_data, 
     +                 lum_mu_data, lum_mub_data,
     +                 b_lum_e, c_lum_e, d_lum_e,
     +                 b_lum_eb, c_lum_eb, d_lum_eb,
     +                 b_lum_mu, c_lum_mu, d_lum_mu,
     +                 b_lum_mub, c_lum_mub, d_lum_mub

      SAVE             last_nu, last_time, last_sn_lum

************************************************************************
*  Initializations
************************************************************************

      DATA             first_call       /.true./

      DATA             t_data           /data_size * 0./,
     +                 lum_e_data       /data_size * 0./,
     +                 lum_eb_data      /data_size * 0./,
     +                 lum_mu_data      /data_size * 0./,
     +                 lum_mub_data     /data_size * 0./
      
      DATA             b_lum_e          /aspline_size * 0./,
     +                 c_lum_e          /aspline_size * 0./,
     +                 d_lum_e          /aspline_size * 0./,
     +                 b_lum_eb         /aspline_size * 0./,
     +                 c_lum_eb         /aspline_size * 0./,
     +                 d_lum_eb         /aspline_size * 0./,
     +                 b_lum_mu         /aspline_size * 0./,
     +                 c_lum_mu         /aspline_size * 0./,
     +                 d_lum_mu         /aspline_size * 0./,
     +                 b_lum_mub        /aspline_size * 0./,
     +                 c_lum_mub        /aspline_size * 0./,
     +                 d_lum_mub        /aspline_size * 0./

      DATA             last_time         /-99999./,
     +                 last_nu           /0/,
     +                 last_sn_lum       /0./

**************************************************************************

      IF (      this_time .EQ. last_time
     +    .AND. nu_type   .EQ. last_nu  ) THEN
         sn_lum = last_sn_lum
         RETURN
      ENDIF

****************************************************************************
****************************************************************************
* Burrows et al. (1992)
****************************************************************************

      IF (sn_model_switch .EQ. 1) THEN

         nu = 0
         sn_lum = 0
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

* functional form of luminosity curves (courtesy Adam Burrows)

         IF (T .LT. TS(NU,2)) THEN
            sn_lum = dsqrt((CO(NU,1) * dexp(T/TAU(NU,1)))**2
     +              + (CO(NU,2) * dexp(T/TAU(NU,2)))**2)
         ELSEIF (T .LT. TS(NU,3)) THEN
            sn_lum = CO(NU,3)
     +              * dexp((T-TS(NU,0))**2*(-0.5d0)/SIG(NU,1)**2)
         ELSEIF ((T .LT. TS(NU,5)) .AND. (nu .NE. 1)) THEN
            sn_lum = 1.d0/(CO(NU,5)*dexp(T/TAU(NU,4)) 
     +              + CO(NU,6) * dexp(T/TAU(NU,5)) * dsin(W(NU)*T 
     +              + THETA(NU)))
     +              + 1.d0 / (CO(NU,4) * dexp(T/TAU(NU,3)))
            sn_lum = 1.d0 / sn_lum
         ELSEIF (T .LT. TS(NU,4)) THEN
            sn_lum = CO(NU,4) * dexp(T/TAU(NU,3))
         ELSEIF (T .LT. TS(NU,5)) THEN
            sn_lum = CO(NU,5) * dexp(T/TAU(NU,4)) 
     +              + CO(NU,6) * dexp(T/TAU(NU,5)) * dsin(W(NU)*T 
     +              + THETA(NU))
         ELSEIF (T .LT. TS(NU,6)) THEN
            sn_lum = CO(NU,7) * dexp((T-TS(NU,6))**2
     +              * (-0.5d0) / SIG(NU,2)**2)
         ELSEIF (T .LT. TS(NU,9)) THEN
            sn_lum = dsqrt((CO(NU,8)*dexp((T-TS(NU,6))**2
     +              * (-0.5d0) / SIG(NU,3)**2))**2 
     +              + (CO(NU,9) * dexp(T/TAU(NU,6)))**2
     +              + (CO(NU,10) / ((1.d0 + T/DELT(NU))**TN(NU)))**2)
         ELSE
            sn_lum = CO(NU,11)
     +              * dexp((T-TS(NU,9))**2 * (-0.5d0)/SIG(NU,4)**2)
         ENDIF

         IF (nu .EQ. 3) sn_lum = sn_lum / 4.d0 ! distinguish "nu_mu" types

         sn_lum = sn_lum * 1.d52   ! x 10^52 ergs/s

* apply luminosity cutoff

         IF (lum_cutoff_flag .EQ. 1) THEN
            IF (this_time .GE. lum_cutoff_time) sn_lum = 0.d0
         ENDIF

*************************************************************************
* Beacom & Vogel model
*************************************************************************

      ELSEIF (sn_model_switch .EQ. 2) THEN

* luminosity: 2-sided exponential (each flavour normalizes to E_bind/6)

         IF (this_time .LE. 0.d0) THEN     ! burst (exponential increase)
            sn_lum = 1.d0 / (tau1 + tau2)
     +              * dexp((this_time - 0.d0) / tau1)
         ELSE                               ! cooling (exponential decay)
            sn_lum = 1.d0 / (tau1 + tau2)
     +              * dexp(-(this_time - 0.d0) / tau2)
         ENDIF

         sn_lum = (E_bind / 6.d0) * sn_lum / ms2s ! ergs/s

* apply luminosity cutoff

         IF (lum_cutoff_flag .EQ. 1) THEN
            IF (this_time .GE. lum_cutoff_time) sn_lum = 0.d0
         ENDIF

*************************************************************************
* Detailed models
*************************************************************************

      ELSEIF (sn_model_switch .GT. 2) THEN

         IF (first_call) THEN

            first_call = .false.
 
            IF (sn_model_switch .EQ. 3) THEN

               lum_data_file = 
     +          'dat/model/Bruenn/15Ms/luminosity.dat'
               lum_data_entries = lum_Bruenn15_entries

            ELSEIF (sn_model_switch .EQ. 4) THEN

               lum_data_file = 
     +          'dat/model/Bruenn/25Ms/luminosity.dat'
               lum_data_entries = lum_Bruenn25_entries

            ELSEIF (sn_model_switch .EQ. 5) THEN

               lum_data_file = 
     +          'dat/model/Mezzacappa/13Ms/luminosity.dat'
               lum_data_entries = lum_Mezz13_entries

            ENDIF

            OPEN (unit = 50, file = lum_data_file, status = 'old')
      
            DO 20 i = 1, lum_data_entries
               READ (50,*) t_data(i), 
     +                     lum_e_data(i), lum_eb_data(i), 
     +                     lum_mu_data(i), lum_mub_data(i)

c               IF (t_data(i) - t_data(i-1) .LE. 0) THEN
c                  write (*,*) 'unordered!: i, t_data(i) =',i, t_data(i)
c               ENDIF

               lum_e_data(i) = dlog10(lum_e_data(i))
               lum_eb_data(i) = dlog10(lum_eb_data(i))
               lum_mu_data(i) = dlog10(lum_mu_data(i))
               lum_mub_data(i) = dlog10(lum_mub_data(i))

20          CONTINUE

            CLOSE (unit = 50)

* nu_e luminosity
            CALL aspline(lum_data_entries,t_data,lum_e_data,
     +                   b_lum_e,c_lum_e,d_lum_e)
 
c            DO i = 1, lum_data_entries
c               write (*,*) t_data(i),lum_e_data(i),
c     +                  b_lum_e(i),c_lum_e(i),d_lum_e(i)
c            END DO

* nu_eb luminosity
            CALL aspline(lum_data_entries,t_data,lum_eb_data,
     +                   b_lum_eb,c_lum_eb,d_lum_eb)

* nu_mu luminosity
            CALL aspline(lum_data_entries,t_data,lum_mu_data,
     +                   b_lum_mu,c_lum_mu,d_lum_mu)

* nu_mub luminosity
            CALL aspline(lum_data_entries,t_data,lum_mub_data,
     +                   b_lum_mub,c_lum_mub,d_lum_mub)

         ENDIF     ! condition for first call

* evaluate the splines for each call to the routine

         IF (nu_type .EQ. id_nu_e) THEN

            CALL aspline_eval(lum_data_entries,
     +                        t_data,lum_e_data,
     +                        b_lum_e,c_lum_e,d_lum_e,
     +                        this_time,lum_fit)

         ELSEIF (nu_type .EQ. id_nu_ebar) THEN

            CALL aspline_eval(lum_data_entries,
     +                        t_data,lum_eb_data,
     +                        b_lum_eb,c_lum_eb,d_lum_eb,
     +                        this_time,lum_fit)

         ELSEIF (     nu_type .EQ. id_nu_mu
     +           .OR. nu_type .EQ. id_nu_tau) THEN

            CALL aspline_eval(lum_data_entries,
     +                        t_data,lum_mu_data,
     +                        b_lum_mu,c_lum_mu,d_lum_mu,
     +                        this_time,lum_fit)

         ELSEIF (     nu_type .EQ. id_nu_mubar
     +           .OR. nu_type .EQ. id_nu_taubar) THEN

            CALL aspline_eval(lum_data_entries,
     +                        t_data,lum_mub_data,
     +                        b_lum_mub,c_lum_mub,d_lum_mub,
     +                        this_time,lum_fit)

         ENDIF   ! check for neutrino type

         IF (      lum_fit .GE. 0.d0
     +       .AND. lum_fit .LT. 1.d60) THEN

            sn_lum = 10.d0**(lum_fit) ! already in ergs/s

         ELSE  ! troubles with Mezzacappa between -90 -> -55 or -38
c            write (*,*) 'bad fit!'
            IF (nu_type .EQ. id_nu_e) THEN
               sn_lum = 5.d50  ! exact # matters here
                               ! bad down to -55
            ELSEIF (nu_type .EQ. id_nu_ebar) THEN
               sn_lum = 8.d44  ! exact # doesn't matter (small contribution)
                               ! bad down to -55
            ELSE
               sn_lum = 2.d44  ! exact # doesn't matter (small contribution)
                               ! bad down to -38
            ENDIF
         ENDIF

c         write (*,*) nu_type, this_time, sn_lum

* apply luminosity cutoff

         IF (lum_cutoff_flag .EQ. 1) THEN
            IF (this_time .GE. lum_cutoff_time) sn_lum = 0.d0
         ENDIF

      ENDIF   ! check for which model

* remember a few things about the last time this routine was called

      last_nu = nu_type
      last_time = this_time
      last_sn_lum = sn_lum

      END
