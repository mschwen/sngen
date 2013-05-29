***********************************************************************
*  Filename: get_flux_index.f
***********************************************************************
*
*
***********************************************************************

      SUBROUTINE get_flux_index(this_energy,this_time,flux_index)

      IMPLICIT none

      INCLUDE 'file_flux_size.inc'

C  Common Block Declarations

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

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

C  Global Variable Declarations

      REAL*8           this_energy, this_time   ! input
      INTEGER          flux_index               ! output

C  Local Variable Declarations

      INTEGER          i                        ! file flux
      INTEGER          this_energy_index, energy_index, time_index
      LOGICAL          found                    ! file flux

      COMMON           /energy_index/
     +                 this_energy_index

*************************************************************************
*  Initializations
*************************************************************************

      energy_index = 0
      time_index   = 0

*************************************************************************
*  Burrows et al. (1992)
*  (time is broken up into burst and cooling regimes)
*************************************************************************

      IF (sn_model_switch .EQ. 1) THEN

* determine energy index
     
         energy_index = nint((this_energy - emin) / estep) *
     +                  time_entries + 1

* determine "this_energy_index" for cross section labels

         this_energy_index = nint((this_energy - emin) / estep) + 1
         
* determine time index

         IF (this_time .LT. time_start2) THEN
            time_index = nint((this_time - time_start1) /
     +                   tstep1)
         ELSE
            time_index = time_entries1
     +                  + nint((this_time - time_start2) / tstep2)
         ENDIF

* determine flux index
              
         flux_index = energy_index + time_index

***************************************************************************
*  Beacom & Vogel (1998)
*  (time and energy intervals are the same for all energy and time)
***************************************************************************

      ELSEIF (sn_model_switch .EQ. 2) THEN

* determine energy index

         energy_index = nint((this_energy - emin) / estep) *
     +                  time_entries + 1

* determine "this_energy_index" for cross section labels

         this_energy_index = nint((this_energy - emin) / estep) + 1

* determine time index

         time_index = nint((this_time - time_start1) / tstep1)

* determine flux index

         flux_index = energy_index + time_index

***************************************************************************
*  Bruenn 25,15 (2000); Mezzacappa (2000)
*  (both time and energy bins are not equal throughout the range)
***************************************************************************

      ELSEIF (sn_model_switch .GE. 3) THEN

* assign energy index for appropriate energy interval

         i = 1
         found = .false.
         DO WHILE (.NOT. found)

            IF (      this_energy .LE. E(i)
     +          .AND. this_energy .GT. E(i-1)) THEN

c               WRITE (*,*) 'this_energy, E(i), E(i-1) =',
c     +                      this_energy, E(i), E(i-1)

               found = .true.
               energy_index = i  ! found the appropriate energy interval
            ELSE
               i = i + 1         ! move on to the next energy interval
            ENDIF

            IF (i .GT. 20 + 1) THEN
c               WRITE (*,*) '*** get_flux_index.f: something went wrong!'
c               WRITE (*,*) '*** get_flux_index.f: energy index, i =', i
            ENDIF
         END DO

         this_energy_index = energy_index

         i = 1
         found = .false.
         DO WHILE (.NOT. found)

c            WRITE (*,*) 'i, this_time, t(i), t(i-1) =',
c     +                   i, this_time, t(i), t(i-1)

            IF (      this_time .LE. t(i)
     +          .AND. this_time .GT. t(i-1)) THEN

c               WRITE (*,*) 'this_time, t(i), t(i-1) =',
c     +                      this_time, t(i), t(i-1)

               found = .true.
               time_index = i  ! found the appropriate energy interval
            ELSE
               i = i + 1       ! move on to the next time interval
            ENDIF

            IF (i .GT. (time_entries + 1)) THEN
c               WRITE (*,*) '*** get_flux_index.f: something went wrong!'
c               WRITE (*,*) '*** get_flux_index.f: time index, i =', i
            ENDIF
         END DO

* how many energy entries do we have to skip over 
* [# time blocks * entries in energy block]

         time_index = (time_index - 1) * file_energy_size 

* determine flux index

         flux_index = energy_index + time_index

      ENDIF

      RETURN
      END
