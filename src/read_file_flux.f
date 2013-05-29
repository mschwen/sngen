*************************************************************************
*  Filename: read_file_flux.f
*************************************************************************
*
*************************************************************************

      SUBROUTINE read_file_flux()

      IMPLICIT none

      INCLUDE 'file_flux_size.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON       /sn_model/
     +             sn_model_switch
      INTEGER      sn_model_switch

C  Local Variable Declarations

      CHARACTER*80 model_file
      INTEGER      i, j, time_step
      REAL*8       t_dum, E_dum, dE_dum
      REAL*8       t(0:file_time_size),
     +             E(0:file_energy_size), dE(file_energy_size),
     +             Le(file_flux_size),
     +             Leb(file_flux_size),
     +             Lmu(file_flux_size)
      COMMON       /t_file/
     +             t
      COMMON       /E_file/
     +             E, dE
      COMMON       /L_file/
     +             Le, Leb, Lmu

      SAVE         /t_file/, /E_file/, /L_file/

*************************************************************************

      WRITE (*,*) '  --> filling flux data array...'

* open reformatted file

      IF (sn_model_switch .EQ. 3) THEN
         model_file = 'dat/model/Bruenn/S15s7b_new.dat'
      ELSEIF (sn_model_switch .EQ. 4) THEN
         model_file = 'dat/model/Bruenn/S25s7b_new.dat'
      ELSEIF (sn_model_switch .EQ. 5) THEN
         model_file = 'dat/model/Mezzacappa/13Ms/Mezzacappa_500.dat'
      ENDIF

      OPEN (unit = 10, file = model_file, status = 'old')

* read in time [s], energy [MeV], d^2/dEdt [1/MeV/s] data (once)

      E(0) = 0.
      t(0) = -200.     ! chose lower than any model start time
      j    = 0

      DO 10 i = 1, file_flux_size

         READ (10,*) t_dum, E_dum, dE_dum, Le(i), Leb(i), Lmu(i)

c         WRITE (*,*) 't_dum, E_dum, dE_dum, Le(i), Leb(i), Lmu(i) =',
c     +                t_dum, E_dum, dE_dum, Le(i), Leb(i), Lmu(i)

* fill energy (E,dE) and time (t) arrays

         time_step = (j * file_energy_size) + 1

         IF (i .EQ. time_step) THEN
            j = j + 1
            t(j) = t_dum * sec2ms     ! convert s to ms in all cases
c            WRITE (*,*) 't =', t(j)
         ENDIF

         IF (i .LE. file_energy_size) THEN
            E(i) = E_dum              ! already in MeV
            dE(i) = dE_dum            ! already in MeV
c            WRITE (*,*) 'E,dE =', E(i), dE(i)
         ENDIF

* sometimes Mezzacappa has negative values that are purely artifacts
* so get rid of them

         IF (Le(i)  .LT. 0.d0) Le(i)  = 0.d0
         IF (Leb(i) .LT. 0.d0) Leb(i) = 0.d0
         IF (Lmu(i) .LT. 0.d0) Lmu(i) = 0.d0

10    CONTINUE  ! loop over all flux entries

      CLOSE (unit = 10)

      RETURN
      END

