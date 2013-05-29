*************************************************************************
*  Filename: count_sum.f
*************************************************************************
*
*  This subroutine will find the neutrino count rate for the chosen
*  interaction(s) over the specified time and energy range. Once that is
*  done, the total number of neutrinos expected in the SNO detector is
*  computed.
*          
*  Inputs: 
*    	o  this interaction id number  (not array)
*  Outputs:
*	o  count total for the specified interaction
*          (then used by sampling routine)
*
*************************************************************************
 
      SUBROUTINE count_sum(this_interaction,count_ave)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'file_flux_size.inc'
      INCLUDE 'convert_units.inc'
      INCLUDE 'time_size.inc'

C  Global Variable Declarations

      INTEGER          this_interaction       ! input
      REAL*8           count_ave              ! output

C  Common Block Declarations

      COMMON           /generator_flag/
     +                 gen_flag, counts_file
      INTEGER          gen_flag
      CHARACTER*40     counts_file

      COMMON           /number_supernova/
     +                 n_supernova
      INTEGER          n_supernova

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /central_stuff/
     +                 central_flag, central_file
      INTEGER          central_flag
      CHARACTER*40     central_file

      COMMON           /energy_limits/ 
     +                 emin, emax, estep,
     +                 dum_emin, dum_emax, dum_estep
      REAL*8           emin, emax, estep
      INTEGER          dum_emin, dum_emax, dum_estep

      COMMON           /time_limits/ 
     +                 time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2,
     +                 dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2
      REAL*8           time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2
      INTEGER          dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2

      COMMON           /sn_time_entries/  
     +                 time_entries1, time_entries2, time_entries
      INTEGER          time_entries1, time_entries2, time_entries

      COMMON           /E_file/
     +                 E, dE
      REAL*8           E(0:file_energy_size), dE(file_energy_size)
      COMMON           /t_file/
     +                 t
      REAL*8           t(0:file_time_size)

* output COMMON blocks

      REAL*8           flux_max, flux_min
      COMMON           /flux_limits/
     +                 flux_max, flux_min

      REAL*8           t_d(time_d_size),
     +                 cr_d(interaction_total,time_d_size),
     +                 c_d(interaction_total,time_d_size),
     +                 cumc_d(interaction_total,time_d_size)
      COMMON           /stats_stuff/
     +                 t_d, cr_d, c_d, cumc_d

C  Local Variable Declarations

      INTEGER          e_dum, t_dum
      REAL*8           time, time_mid, time_delay, time_d, time_d_stop
      REAL*8           energy
      REAL*8           sn_flux, det_flux,   ! function
     +                 this_sn_flux, this_det_flux

      REAL*8           central,             ! function
     +                 this_central, central_max, energy_max

      INTEGER          time_index1, time_index2, time_index
      INTEGER          time_d_index1, time_d_index2, time_d_index3,
     +                 time_d_index, time_d_index_max

      REAL*8           countrate, d_countrate,
     +                 tstep, factor

      REAL*8           this_countrate, last_countrate
      REAL*8           this_time, last_time
      REAL*8           countrate_ave, count, count_tot
      REAL*8           time_interval

      LOGICAL          first_call
      SAVE             first_call

      INTEGER          i, j, k

      DATA             first_call /.true./     ! set true to start

**************************************************************************

* check to see if we need to sum counts 
* (sampling may already be done and we just need to output SNOMAN files)

      IF (gen_flag .EQ. 3) THEN

         OPEN (unit=25, file=counts_file, status='unknown')

* read counts file until the line with the specified interaction
* (count_ave should then have the correct counts for this interaction)

         WRITE (*,*)
         WRITE (*,*) '  Reading supernova neutrino counts for interactio
     +n ', this_interaction,'...'

         DO i = 1, this_interaction
            READ (25,*) j, count_ave    ! interaction id, counts
c            WRITE (*,*) 'count_ave = ', count_ave
         END DO
         CLOSE (unit = 25)

         count_ave = count_ave * dble(n_supernova)

         WRITE (*,5) this_interaction, count_ave
5        FORMAT (3x,'For interaction ',i2,' the total number of counts i
     +s ' ,f9.3)
         RETURN
      ENDIF

**************************************************************************
*  Initializations
**************************************************************************

      j = 0
      time_index1 = 0
      time_index2 = 0
      time_index  = 0
      time_d_index = 0
      time_d_index1 = 0
      time_d_index2 = 0
      time_d_index3 = 0
      time_d_index_max = 0

      count_tot = 0.d0

      DO i = 1, time_d_size
         t_d(i) = 0.d0
      END DO

      IF (first_call) THEN
         DO i = 1, interaction_total
            DO k = 1, time_d_size
               cr_d(i,k) = 0.d0
               c_d(i,k) = 0.d0
               cumc_d(i,k) = 0.d0
            END DO
         END DO
         IF (sn_model_switch .LE. 2) first_call = .false.
      ENDIF

      flux_min = 1000.
      flux_max = -999.                                                       

      last_countrate = -1.                       

**************************************************************************
*  Start the integration over time and energy
**************************************************************************     

      WRITE (*,*)
      WRITE (*,*) '  Calculating supernova neutrino counts for interacti
     +on ', this_interaction,'...'

**************************************************************************
**************************************************************************
*  Burrows et al. (1992)
**************************************************************************

      IF (sn_model_switch .EQ. 1) THEN  ! Burrows

**************************************************************************
*  Calculate the count rate for the various interaction(s), one
*  interaction at a time.
*
*  Start with det_flux in units of [counts/(MeV*s)]:
*       o  integrate over energy to get the count rate as a function
*          of time
*       o  integrate over time to get the total number of counts
*
*  Note:  last_countrate < 0 allows the loop to skip over the first
*         calculation of the total count (needs cr(i) + cr(i-1))
**************************************************************************

**************************************************************************
*  Burst:
*  Integrate "counts/(MeV*s)" over energy to get count rate as
*  function of time
**************************************************************************
* loop over supernova neutrino time (not delayed)

         DO 100 t_dum = dum_tmin1, dum_tmax1, dum_tstep1
            time = dble(t_dum) * tstep1
            time_mid = dble(t_dum - 0.5d0) * tstep1  ! (this - last)

* time index for count rate (undelayed)

            time_index1 = nint((time - time_start1) / tstep1 + 1)
            time_index = time_index1

* initialize count rate counter for each time

            countrate = 0.d0

* loop over supernova neutrino energy

            DO 110 e_dum = dum_emin, dum_emax, dum_estep
               energy = dble(e_dum) * estep

* get supernova source and detector flux values

               this_sn_flux = sn_flux(this_interaction,energy,time)
               this_det_flux = det_flux(this_interaction,energy,
     +                                  this_sn_flux)

* look for maximum in burst stage (might as well look for minimum, too,
* but it will most likely be found in the cooling phase)
* [output using COMMON blocks and used in defining flux sampling bounds]

               IF (this_det_flux .GT. flux_max)
     +          flux_max = this_det_flux
               IF (this_det_flux .LT. flux_min) ! prob won't find here 
     +          flux_min = this_det_flux

* differential count rate for this energy and this time;
* also account for the number of supernova bursts

               d_countrate = this_det_flux * estep * dble(n_supernova)

* calculate energy-dependent time delay for massive neutrinos
* (if neutrino masses are zero, function just returns undelayed time)

               time_d = time_delay(this_interaction,energy,time)

               IF (time_d .EQ. time) THEN  ! not delayed

* sum differential count rate (ie., integrate over energy)

                  countrate = countrate + d_countrate

* fill count rate arrays for delayed time and label accordingly

               ELSEIF (time_d .NE. time) THEN

c                  WRITE (*,*) 'We''re delayed: time,time_d =',
c     +                         time, time_d

c                   IF (      time_d .GE. 1000.d0
c     +                 .and. time_d .LT. 1020.d0)
c     +              write (*,*) 'time,energy,time_d=',
c     +                           time,energy,time_d

* a few low energy events will be drastically delayed; ignore them for now

                  time_d_stop = 5.d0 * time_stop2

                  IF (time_d .GT. time_d_stop) GOTO 110

* determine the time index of the delayed time
* (need to consider both burst and cooling regimes since neutrinos
*  produced in the burst phase may be shifted to the cooling phase)

                  IF (time_d .GT. time_stop1) THEN

* need to account for the different integration steps used in the burst and
* cooling stages (hence, the "factor")
                     factor = tstep1 / tstep2 
                     time_d_index1 = nint(time_stop1 - time_start1) 
     +                              / tstep1 + 1
                     IF (time_d .GT. time_stop2) THEN
                        time_d_index2 = nint(time_stop2 - time_start2)
     +                                / tstep2 + 1
                        time_d_index3 = nint(time_d - time_stop2) 
     +                                / tstep2  ! no "+1" becuase not using start
                     ELSEIF (time_d .LE. time_stop2) THEN
                        time_d_index2 = nint(time_d - time_start2)
     +                                / tstep2 + 1
                        time_d_index3 = 0
                     ENDIF
                  ELSEIF (time_d .LE. time_stop1) THEN
                     factor = 1.d0
                     time_d_index1 = nint(time_d - time_start1)
     +                             / tstep1 + 1
                     time_d_index2 = 0
                     time_d_index3 = 0
                  ENDIF

                  time_d_index =  time_d_index1
     +                          + time_d_index2 + time_d_index3

* quantize the delayed time in terms of the various time steps

                  t_d(time_d_index) =  time_start1 
     +                               + (time_d_index1 - 1) * tstep1
     +                               +  time_d_index2 * tstep2
     +                               +  time_d_index3 * tstep2

* we want the mid-point of this time bin (undelayed already calculated)

                  IF (time_d_index2 .EQ. 0) THEN      ! burst phase
                     tstep = tstep1
                  ELSEIF (time_d_index2 .GT. 0) THEN  ! cooling phase
                     tstep = tstep2
                  ENDIF

                  t_d(time_d_index) =  t_d(time_d_index)  ! (this - last)
     +                               - 1.5d0 * tstep

* assign the differential count rates at a particular neutrino energy to
* the appropriate time bins

                  cr_d(this_interaction,time_d_index) = 
     +              cr_d(this_interaction,time_d_index)
     +             + d_countrate * factor

c                   IF (      time_d .GE. 1000.d0
c     +                 .and. time_d .LE. 1020.d0)
c     +              write (*,*) 'time,energy,time_d,time_d_index,
c     + t_d(time_d_index),factor,d_countrate,
c     + cr_d(this_interaction,time_d_index)=',time,energy,time_d,
c     + time_d_index,t_d(time_d_index),factor,d_countrate,
c     + cr_d(this_interaction,time_d_index)

* write debugging info at certain time (values will change because time delay
* depends on energy)

c                  IF (time .EQ. 500.) THEN
c                     write (*,*) 'time,time_d,time_d_index,t_d,cr_d =',
c     +                            time,time_d,time_d_index,
c     +                            t_d(time_d_index),
c     +                            cr_d(this_interaction,time_d_index)
c                  ENDIF

               ENDIF  ! check whether delayed time or not
110         CONTINUE  ! energy loop

**************************************************************************
*  Burst:  [massless neutrinos]
*  Integrate count rate for all times (average count rate over time bin)
*  [Remember to convert time to seconds]
*
* if we're not dealing with massive neutrinos, then save time and fill
* count rate arrays the easy way without having to bin (compares the
* delayed and undelayed times at the highest energy to establish whether
* we have massive neutrinos or not)
**************************************************************************

            IF (time_d .EQ. time) THEN
               t_d(time_index) = time_mid
               cr_d(this_interaction,time_index) = countrate

               this_countrate = countrate
               this_time = time

               IF (time_index .GT. 1) THEN  ! skip over first time bin
                  countrate_ave = 0.5d0 * (this_countrate + 
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s) 
                  count_tot = count_tot + count 

                  c_d(this_interaction,time_index) = count  ! raw counts
                  cumc_d(this_interaction,time_index) = count_tot  ! cumulative
               ENDIF

               last_countrate = this_countrate
               last_time = this_time
            ENDIF  ! time integration for the no time delay case
100      CONTINUE  ! (undelayed) burst time loop

**************************************************************************
*  Burst: [massive neutrinos]
*  Integrate count rate for all times (average count rate over time bin)
*  (loop over supernova neutrino time (not delayed))
*  [Remember to convert time to seconds]
**************************************************************************

         IF (time_d .NE. time) THEN

            DO 120 t_dum = dum_tmin1, dum_tmax1, dum_tstep1
               time = dble(t_dum) * tstep1

* time index for count rate
* (want to compare sequential indices, so use undelayed time)

               time_index1 = nint((time - time_start1) / tstep1 + 1)
               time_index = time_index1

               this_countrate = cr_d(this_interaction,time_index)
               this_time = time

               IF (time_index .GT. 1) THEN  ! skip over first time bin
                  countrate_ave = 0.5d0 * (this_countrate +
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s)
                  count_tot = count_tot + count

c                  write (*,*) 'time =', time
c                  write (*,*) 'count, count_tot=', count, count_tot

                  c_d(this_interaction,time_index) = count  ! raw counts
                  cumc_d(this_interaction,time_index) = count_tot  ! cumulative
               ENDIF  ! time integration

               last_countrate = this_countrate
               last_time = this_time
120         CONTINUE ! (undelayed) burst time loop
         ENDIF       ! check for delayed time (massive neutrinos)

**************************************************************************
*  Cooling:
*  Integrate d2N/dEdt [counts/(MeV*s)] over energy to get count rate as
*  function of time
**************************************************************************

* loop over supernova neutrino time (not delayed)

         DO 150 t_dum = dum_tmin2, dum_tmax2, dum_tstep2
            time = dble(t_dum) * tstep2
            time_mid = dble(t_dum - 0.5d0) * tstep2  ! (this - last)

* set up time index labels for group count rate calculation

            time_index2 = nint((time - time_start2) / tstep2 + 1)
            time_index = time_index1 + time_index2

* initialize count rate counter

            countrate = 0.d0

* loop over supernova neutrino energy

            DO 160 e_dum = dum_emin, dum_emax, dum_estep
               energy = dble(e_dum) * estep

* get detector flux for given time and energy

               this_sn_flux = sn_flux(this_interaction,energy,time)
               this_det_flux = det_flux(this_interaction,energy,
     +                                  this_sn_flux)

* look for minimum in cooling stage (might as well look for maximum, too)

               IF (this_det_flux .GT. flux_max)  ! prob won't find here
     +          flux_max = this_det_flux
               IF (this_det_flux .LT. flux_min)
     +          flux_min = this_det_flux

* differential count rate

               d_countrate = this_det_flux * estep * dble(n_supernova)

* calculate energy-dependent time delay for massive neutrinos
* (if neutrino masses are zero, function just returns undelayed time)

               time_d = time_delay(this_interaction,energy,time)

               IF (time_d .EQ. time) THEN

* sum differential count rate (ie., integrate over energy)

                  countrate = countrate + d_countrate

* fill count rate arrays for delayed time and label accordingly

               ELSEIF (time_d .NE. time) THEN

c                  WRITE (*,*) 'We''re delayed: time,time_d =',
c     +                         time, time_d

                  IF (time_d .GT. time_d_stop) GOTO 160

* determine the time index of the delayed time
* (no "factor" required since all events here will be in cooling phase)

                  IF (time_d .GT. time_stop1) THEN
                     time_d_index1 = nint(time_stop1 - time_start1)
     +                              / tstep1 + 1
                     IF (time_d .GT. time_stop2) THEN
                        time_d_index2 = nint(time_stop2 - time_start2)
     +                              / tstep2 + 1
                        time_d_index3 = nint(time_d - time_stop2)
     +                              / tstep2  ! no "+1" because not using start
                     ELSEIF (time_d .LE. time_stop2) THEN
                        time_d_index2 = nint(time_d - time_start2)
     +                              / tstep2 + 1
                        time_d_index3 = 0
                     ENDIF
                  ELSEIF (time_d .LE. time_stop1) THEN
                     time_d_index1 = nint(time_d - time_start1)
     +                           / tstep1 + 1
                     time_d_index2 = 0
                     time_d_index3 = 0
                  ENDIF

                  time_d_index =  time_d_index1
     +                          + time_d_index2 + time_d_index3

* quantize the delayed time in terms of the various time steps

                  t_d(time_d_index) =  time_start1
     +                               + (time_d_index1 - 1) * tstep1
     +                               +  time_d_index2 * tstep2
     +                               +  time_d_index3 * tstep2

* we want the mid-point of this time bin (undelayed already calculated)

                  IF (time_d_index2 .EQ. 0) THEN
                     tstep = tstep1
                  ELSEIF (time_d_index2 .GT. 0) THEN
                     tstep = tstep2
                  ENDIF

                  t_d(time_d_index) =  t_d(time_d_index)  ! (this - last)
     +                               - 0.5d0 * tstep 

                  cr_d(this_interaction,time_d_index) =
     +              cr_d(this_interaction,time_d_index)
     +             + d_countrate

c                   IF (      time_d .GE. 1000.d0
c     +                 .and. time_d .LT. 1040.d0)
c     +              write (*,*) 'time,energy,time_d,time_d_index,
c     + t_d(time_d_index),factor,d_countrate,
c     + cr_d(this_interaction,time_d_index)=',time,energy,time_d,
c     + time_d_index,t_d(time_d_index),factor,d_countrate,
c     + cr_d(this_interaction,time_d_index)

               ENDIF  ! check whether time is delayed or not
160         CONTINUE  ! energy loop

**************************************************************************
*  Cooling:  [massless neutrinos]
*  Integrate count rate, dN/dt, over all times
*  (average count rate over time bin)
*  [Remember to convert time to seconds]
*
* if we're not dealing with massive neutrinos, then save time and fill
* count rate arrays the easy way without having to bin (compares the
* delayed and undelayed times at the highest energy to establish whether
* we have massive neutrinos or not)
**************************************************************************

            IF (time_d .EQ. time) THEN
               t_d(time_index) = time_mid
               cr_d(this_interaction,time_index) = countrate

               this_countrate = countrate
               this_time = time

               IF (time_index .GT. 1) THEN
                  countrate_ave = 0.5d0 * (this_countrate + 
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s)
                  count_tot = count_tot + count

                  IF (time_d .EQ. time) THEN
                     c_d(this_interaction,time_index) = count  ! raw counts
                     cumc_d(this_interaction,time_index) = count_tot  ! cumulative
                  ENDIF
               ENDIF

               last_countrate = this_countrate
               last_time = this_time
            ENDIF  ! time integration for undelayed time case
150      CONTINUE  ! (undelayed) cooling time loop

**************************************************************************
*  Cooling: [massive neutrinos]
*  Integrate count rate for all times (average count rate over time bin)
*  (loop over supernova neutrino time (not delayed))
*  [Remember to convert time to seconds]
**************************************************************************

         IF (time_d .NE. time) THEN

            DO 170 t_dum = dum_tmin2, dum_tmax2, dum_tstep2
               time = dble(t_dum) * tstep2

* time index for count rate
* (want to compare sequential indices, so use undelayed time)

               time_index2 = nint((time - time_start2) / tstep2 + 1)
               time_index = time_index1 + time_index2

               this_countrate = cr_d(this_interaction,time_index)
               this_time = time

               IF (time_index .GT. 1) THEN  ! skip over first time bin
                  countrate_ave = 0.5d0 * (this_countrate +
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s)
                  count_tot = count_tot + count

c                  write (*,*) 'time =', time
c                  write (*,*) 'count, count_tot=', count, count_tot

                  c_d(this_interaction,time_index) = count  ! raw counts
                  cumc_d(this_interaction,time_index) = count_tot  ! cumulative
               ENDIF

               last_countrate = this_countrate
               last_time = this_time
170         CONTINUE ! (undelayed) cooling time loop
         ENDIF       ! check for delayed time (massive neutrinos)

         count_ave = count_tot

         WRITE (*,180) this_interaction, count_ave
180      FORMAT (3x,'For interaction ',i2,' the total number of counts i
     +s ' ,f9.3)

c         write (*,*) 'flux_max = ', flux_max  ! shows max doesn't change
c         write (*,*) 'flux_min = ', flux_min  ! shows min is found here

*************************************************************************
*************************************************************************
*  Beacom & Vogel (1998)
*************************************************************************

      ELSEIF (sn_model_switch .EQ. 2) THEN ! Beacom & Vogel

* loop over supernova neutrino time (not delayed)

         DO 200 t_dum = dum_tmin1, dum_tmax1, dum_tstep1
            time = dble(t_dum) * tstep1
            time_mid = dble(t_dum - 0.5d0) * tstep1  ! (this - last)

* time index for count rate (undelayed)

            time_index = nint((time - time_start1) / tstep1 + 1)

* initialize count rate counter for each time

            countrate = 0.d0

* loop over supernova neutrino energy

            DO 210 e_dum = dum_emin, dum_emax, dum_estep
               energy = dble(e_dum) * estep

* get supernova source and detector flux values

               this_sn_flux = sn_flux(this_interaction,energy,time)
               this_det_flux = det_flux(this_interaction,energy,
     +                                  this_sn_flux)

* look for maximum in burst stage (might as well look for minimum, too)

               IF (this_det_flux .GT. flux_max)
     +          flux_max = this_det_flux
               IF (this_det_flux .LT. flux_min) ! prob won't find here 
     +          flux_min = this_det_flux

* differential count rate

               d_countrate = this_det_flux * estep * dble(n_supernova)

* calculate energy-dependent time delay for massive neutrinos
* (if neutrino masses are zero, function just returns undelayed time)

               time_d = time_delay(this_interaction,energy,time)

               IF (time_d .EQ. time) THEN  ! not delayed
 
* sum differential count rate (ie., integrate over energy)

                  countrate = countrate + d_countrate

* fill count rate arrays for delayed time and label accordingly

               ELSEIF (time_d .NE. time) THEN

c                  WRITE (*,*) 'We''re delayed: time,time_d =',
c     +                         time, time_d

* a few low energy events will be drastically delayed; ignore them for now

                  time_d_stop = 5.d0 * time_stop1

                  IF (time_d .GT. time_d_stop) GOTO 210

* determine the time index of the delayed time

                  IF (time_d .GT. time_stop1) THEN
                     time_d_index1 = nint(time_stop1 - time_start1) 
     +                              / tstep1 + 1
                     time_d_index2 = nint(time_d - time_stop1) 
     +                              / tstep1  ! no "+1" becuase not using start
                  ELSEIF (time_d .LE. time_stop1) THEN
                     time_d_index1 = nint(time_d - time_start1)
     +                             / tstep1 + 1
                     time_d_index2 = 0
                  ENDIF

                  time_d_index =  time_d_index1
     +                          + time_d_index2 + time_d_index3

                  IF (time_d_index .GT. time_d_index_max) THEN
                     time_d_index_max = time_d_index
c                     write (*,*) 'time_d_index_max,time,energy =',
c     +               time_d_index_max,time,energy
                  ENDIF
 
* quantize the delayed time in terms of the various time steps

                  t_d(time_d_index) =  time_start1 
     +                               + (time_d_index1 - 1) * tstep1
     +                               +  time_d_index2 * tstep1

* we want the mid-point of this time bin (undelayed already calculated)

                  t_d(time_d_index) =  t_d(time_d_index)  ! (this - last)
     +                               - 0.5d0 * tstep1

* assign the differential count rates at a particular neutrino energy to
* the appropriate time bins

                  cr_d(this_interaction,time_d_index) = 
     +              cr_d(this_interaction,time_d_index)
     +             + d_countrate

* write debugging info at certain time (values will change because time delay
* depends on energy)

c                  IF (time .EQ. time_stop1) THEN
c                     write (*,*) 'time,time_d,time_d_index,time_d_index_
c     +max,t_d,cr_d =',
c     +                            time,time_d,time_d_index,
c     +                            time_d_index_max,
c     +                            t_d(time_d_index),
c     +                            cr_d(this_interaction,time_d_index)
c                  ENDIF
               ENDIF  ! delayed time stuff
210         CONTINUE  ! energy loop

**************************************************************************
*  [Massless neutrinos]
*  Integrate count rate for all times (average count rate over time bin)
*  [Remember to convert time to seconds]
*
* If we're not dealing with massive neutrinos, then save time and fill
* count rate arrays the easy way without having to bin (compares the
* delayed and undelayed times at the highest energy to establish whether
* we have massive neutrinos or not)
**************************************************************************

            IF (time_d .EQ. time) THEN
               t_d(time_index) = time_mid
               cr_d(this_interaction,time_index) = countrate

               this_countrate = countrate
               this_time = time

               IF (time_index .GT. 1) THEN     ! skip first time bin
                  countrate_ave = 0.5d0 * (this_countrate + 
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s) 
                  count_tot = count_tot + count 

                  c_d(this_interaction,time_index) = count  ! raw counts
                  cumc_d(this_interaction,time_index) = count_tot  ! cumulative
               ENDIF

               last_countrate = this_countrate
               last_time = this_time
            ENDIF   ! time integration for the no time delay case
200      CONTINUE   ! (undelayed) time loop

**************************************************************************
*  [Massive neutrinos]
*  Integrate count rate for all times (average count rate over time bin)
*  (loop over supernova neutrino time (not delayed))
*  [Remember to convert time to seconds]
**************************************************************************

* time integration for the time delayed case (loop over supernova neutrino 
* time (not delayed))

         IF (time_d .NE. time) THEN

c            write (*,*) 'time_d_index_max =', time_d_index_max

            DO 220 t_dum = dum_tmin1, time_d_index_max, dum_tstep1
               time = dble(t_dum) * tstep1

* time index for count rate (undelayed)

               time_index = nint((time - time_start1) / tstep1 + 1)

               this_countrate = cr_d(this_interaction,time_index)
               this_time = time

               IF (time_index .GT. 1) THEN  ! skip over first time bin
                  countrate_ave = 0.5d0 * (this_countrate +
     +                                     last_countrate)
                  time_interval = (this_time - last_time)

                  count = countrate_ave * (time_interval * ms2s)
                  count_tot = count_tot + count

c                  write (*,*) 'time =', time
c                  write (*,*) 'count, count_tot=', count, count_tot

                  c_d(this_interaction,time_index) = count  ! raw counts
                  cumc_d(this_interaction,time_index) = count_tot  ! cumulative
               ENDIF

               last_countrate = this_countrate
               last_time = this_time
220         CONTINUE ! (undelayed) time loop
         ENDIF       ! delayed time

         count_ave = count_tot

         WRITE (*,215) this_interaction, count_ave
215      FORMAT (3x,'For interaction ',i2,' the total number of counts i
     +s ' ,f9.3)

c         write (*,*) 'flux_max = ', flux_max  ! shows max doesn't change
c         write (*,*) 'flux_min = ', flux_min  ! shows min is found here

*************************************************************************
*************************************************************************
*  Breunn 25 Ms, 15Ms (2000); Mezzacappa 13Ms (2000)
*************************************************************************

      ELSEIF (     sn_model_switch .EQ. 3   ! Bruenn_15
     +        .OR. sn_model_switch .EQ. 4   ! Bruenn_25
     +        .OR. sn_model_switch .EQ. 5)  ! Mezzacappa
     +   THEN

         IF (first_call) THEN
            CALL read_file_flux()
            first_call = .false.
         ENDIF

         DO 300 t_dum = 1, file_time_size

* initialize count rate counter for each time

            countrate = 0.d0

            time = t(t_dum)

            DO 310 E_dum = 1, file_energy_size

               energy = E(E_dum)

* get supernova and detector flux values

c               WRITE (*,*) 't_dum, time =', t_dum, time
c               WRITE (*,*) 'E_dum, energy =', E_dum, energy

               this_sn_flux = sn_flux(this_interaction,energy,time)
               this_det_flux = det_flux(this_interaction,energy,
     +                                  this_sn_flux)

* look for both maximum and minimum over all times 

               IF (this_det_flux .GT. flux_max)
     +          flux_max = this_det_flux
               IF (this_det_flux .LT. flux_min)
     +          flux_min = this_det_flux

c               WRITE (*,*) 'sn_flux =', this_sn_flux
c               WRITE (*,*) 'det_flux =', this_det_flux


* differential count rate

               d_countrate = this_det_flux * dE(E_dum)
     +                      * dble(n_supernova)

c               WRITE (*,*) 'dE(E_dum),d_countrate = ',
c     +                      dE(E_dum),d_countrate

* sum differential count rate (ie., perform integration over energy)

               countrate = countrate + d_countrate

c               WRITE (*,*) 'countrate =', countrate

* delay time for count rate with massive neutrinos
* (forget for now)

               time_d = time

310         CONTINUE  ! energy loop

* if we're not dealing with massive neutrinos, then save time and fill
* count rate arrays the easy way without having to bin (compares the
* delayed and undelayed times at the highest energy to establish whether
* we have massive neutrinos or not)

            IF (time_d .EQ. time) THEN
               t_d(t_dum) = time
               cr_d(this_interaction,t_dum) = countrate
            ENDIF

**************************************************************************

            this_countrate = countrate
            this_time = time

            IF (last_countrate .GE. 0.) THEN
               countrate_ave = 0.5d0 * (this_countrate +
     +                                  last_countrate)
               time_interval = (this_time - last_time)

c               WRITE (*,*) 'time_interval =', time_interval

               count = countrate_ave * (time_interval * ms2s)

* sum differential counts (ie., perform integration over time)

               count_tot = count_tot + count

               IF (time_d .EQ. time) THEN
                  c_d(this_interaction,t_dum) = count  ! raw counts
                  cumc_d(this_interaction,t_dum) = count_tot  ! cumulative
               ENDIF
            ENDIF

            last_countrate = this_countrate
            last_time = this_time

300      CONTINUE    ! (undelayed) time loop

* multiply final interaction count total by the number of supernovae

         count_ave = count_tot

         WRITE (*,315) this_interaction, count_ave
315      FORMAT (3x,'For interaction ',i2,' the total number of counts i
     +s ' ,f9.3)

      ENDIF  ! SN models

*************************************************************************

* consider central energies at end of routine since it screws up accounting
* otherwise (eg., xsect routines are expecting one entire pass through all
* energies to fill arrays for faster access)

      IF (central_flag .EQ. 1) THEN

         central_max = -999.
         energy_max = 0.

* loop over supernova neutrino time (not delayed)

         DO 500 t_dum = dum_tmin1, dum_tmax1, dum_tstep1
            time = dble(t_dum) * tstep1

            DO 510 e_dum = dum_emin, dum_emax, dum_estep
               energy = dble(e_dum) * estep

               this_central = central(this_interaction,energy,time)
               IF (this_central .GT. central_max) THEN
                  central_max = this_central
                  energy_max = energy
               ENDIF

510         CONTINUE
500      CONTINUE

         WRITE (75,*) this_interaction, energy_max, central_max
      ENDIF  ! check for central flag
      
      RETURN
      END
