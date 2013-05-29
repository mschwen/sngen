***********************************************************************
*  Filename: sample_nu.f
***********************************************************************
*
*  This routine will sample the normalized detected neutrino
*  as calculated from the previous integration routines. 
*  The energy is converted directly to the particle detected.
*  The method of sampling is the hit and miss rejection method
*  The random number generator is "ran3" from Numerical Recipes.       
*
*  Inputs:
*	o  this interaction id number (not array)
*	o  number of sample events to produce
*  Outputs:
*       o  sampled neutrino time, energy (in "sample_file")
*
*********************************************************************** 

      SUBROUTINE sample_nu(interaction,nevent)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'sample_count_size.inc'
      INCLUDE 'file_flux_size.inc'

C  Common Block Declarations

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

      COMMON           /t_file/
     +                 t
      REAL*8           t(0:file_time_size)

      COMMON           /flux_limits/
     +                 flux_max, flux_min
      REAL*8           flux_max, flux_min

      COMMON           /nu_sample_file/
     +                 nu_sample_file
      CHARACTER*40     nu_sample_file

      COMMON           /sample_efficiency/
     +                 efficiency_flag
      INTEGER          efficiency_flag

* output COMMON blocks

      REAL*8           sample_E_nu(interaction_total,sample_count_size)

      COMMON           /sample_stuff/
     +                 sample_E_nu

C  Global Variable Declarations

      INTEGER          interaction, nevent     ! input

C  Local Variable Declarations

      LOGICAL          first_call
      INTEGER          i, j
      INTEGER          event_id

      DOUBLE PRECISION E_thres

      REAL*8           time_min, time_max, time_ran

      REAL*8           sn_flux, det_flux,    ! function
     +                 this_sn_flux,
     +                 flux_ran, flux_test

      DOUBLE PRECISION ran_val
      REAL*8           energy_ran

      REAL*8           time_delay
      REAL*8           sample_time_d(sample_count_size), 
     +                 sample_energy(sample_count_size)

      INTEGER          count, event_count

      INTEGER          ncall
      REAL*8           efficiency

      REAL*8           first_time_id, last_time_id
      INTEGER          ptr
      REAL*8           time_tmp, energy_tmp      

      SAVE             first_call, count
      DATA             first_call /.true./
      DATA             count /0/

**********************************************************************
*  Initializations
**********************************************************************

      ncall = 0
      event_id = 0

      event_count = 0

      DO 12 i = 1, sample_count_size
         sample_time_d(i) = 0.
         sample_energy(i) = 0.
12    CONTINUE

      IF (first_call) THEN
         DO i = 1, interaction_total
            DO j = 1, sample_count_size
               sample_E_nu(i,j) = 0.d0
            END DO
         END DO
         first_call = .false.
      ENDIF

*************************************************************************

      WRITE (*,5) interaction, nevent
5     FORMAT (1x,'  Sampling interaction',i3,' to produce SN neutrinos f
     +or',i7,' events')

************************************************************************
*  Begin sampling the flux distribution by sampling energy to define
*  which flux vs. time distribution to further sample from.
*  Will stop when the requested number of sample events ("nevents") is reached
************************************************************************

      CALL file_open(91,nu_sample_file,interaction)

      DO WHILE (event_id .LT. nevent)

         ncall = ncall + 1

c         write (*,*) 'ncall = ', ncall

************************************************************************
*  Energy sampling.  "energy_ran" will give a uniform sampling over the
*  entire energy interval
************************************************************************

* get interaction energy threshold

         CALL get_threshold(interaction,E_thres)

* sample energy
* (make sure we don't round down to below E_thres)

         energy_ran = E_thres   ! initial value for comparison loop

         IF (emin .LT. E_thres) THEN
            DOWHILE ((energy_ran - E_thres) .LT. estep)
               energy_ran = ran_val(E_thres,emax)
            END DO

* make sure we don't round down to below emin

         ELSE
            DOWHILE ((energy_ran - emin) .LT. estep)
               energy_ran = ran_val(emin,emax)
            END DO
         ENDIF

c         WRITE (*,*) 'energy_ran =', energy_ran

* just to make sure...

         IF (energy_ran .LT. E_thres) THEN
            WRITE (*,*) '*** sample_nu.f: Something bad happened!'
            WRITE (*,*) '*** sample_nu.f: --> energy_ran < E_thres:',
     +                   energy_ran
         ENDIF

***********************************************************************
*  Time sampling. "time_ran" will give a uniform sampling over the full
*  time delayed interval. 
***********************************************************************

* assign maximum and minimum time values (do not delay because the
* sn_flux and det_flux routines are not expecting delayed times)

         time_min = time_start1
         time_max = time_stop2

* find a random time between the minimum and maximum allowed
* (undelayed) time values

         time_ran = ran_val(time_min,time_max)

c         WRITE (*,*) 'time_ran =', time_ran

***********************************************************************
*  Flux sampling. "flux_ran" will give a uniform sampling from the
*  minimum flux to the maximum flux in the file. 
*  Note: max,min flux limits for all time and energy determined in 
*  "count_sum" routine; in order to sample fairly, these limits need to be
*  the same for all trials.
***********************************************************************

         flux_ran = ran_val(flux_min,flux_max)

***********************************************************************
*  Reads off the value of the flux at the energy and (undelayed) time 
*  and tests to see if it is less than the sample value. Reject if
*  greater, otherwise accept.
***********************************************************************

         this_sn_flux = sn_flux(interaction,energy_ran,time_ran)
         flux_test = det_flux(interaction,energy_ran,this_sn_flux)

         IF (flux_ran .LE. flux_test) THEN  

c            write (*,*) 'found one!'
            event_id = event_id + 1
            count = count + 1
            event_count = count

* fill energy and (delayed) time arrays that can be sorted

            sample_time_d(event_id) = time_delay(interaction,
     +                                energy_ran,time_ran)

            IF (energy_ran .GE. 170.d0) energy_ran = 169.d0 ! xsect limit
            sample_energy(event_id) = energy_ran

         ENDIF  ! check for successful sample
      END DO    ! do while need more sampled particles

*************************************************************************
*  Sort sampled events according to time using a "selection sort"
*  algorithm (this helps SNOMAN merge the reaction events into a single
*  burst file)
*************************************************************************

      last_time_id = event_id
      DO 40 i = 1, last_time_id-1
         ptr = i
         first_time_id = i + 1
         DO 41 j = first_time_id, last_time_id
            IF (sample_time_d(j) .LT. sample_time_d(ptr)) ptr = j
41       CONTINUE
         time_tmp = sample_time_d(i)
         energy_tmp = sample_energy(i)

         sample_time_d(i) = sample_time_d(ptr)
         sample_energy(i) = sample_energy(ptr)

         sample_time_d(ptr) = time_tmp
         sample_energy(ptr) = energy_tmp

         sample_E_nu(interaction,i) = sample_energy(i)

         WRITE (91,45) sample_time_d(i), sample_energy(i)

40    CONTINUE

      sample_E_nu(interaction,last_time_id) = 
     + sample_energy(last_time_id)

      WRITE (91,45) sample_time_d(last_time_id), 
     +              sample_energy(last_time_id)
45    FORMAT (f15.8,f14.8)

      CLOSE (unit = 91)

************************************************************************
*  Calculate sampling efficiency
************************************************************************

      IF (efficiency_flag .EQ. 1) THEN
         efficiency = float(nevent) / float(ncall)

         WRITE (*,55) efficiency
55       FORMAT (1x,'The sampling efficiency is ',e12.5)
      ENDIF

************************************************************************

      RETURN
      END
