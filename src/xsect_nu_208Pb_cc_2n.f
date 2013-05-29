***********************************************************************
*  Filename: xsect_nu_208Pb_cc_2n.f
***********************************************************************
*
*  Inputs:
*	o  interaction id number (not array)
*	o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy
*
***********************************************************************
     
      SUBROUTINE xsect_nu_208Pb_cc_2n(interaction,E_nu,xsect_total)

      IMPLICIT none

      INCLUDE 'energy_size.inc'
      INCLUDE 'xsect_size.inc'
      INCLUDE 'aspline_size.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON           /generator_flag/
     +                 gen_flag, counts_file
      INTEGER          gen_flag
      CHARACTER*40     counts_file

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

C  Global Variable Declarations

      INTEGER interaction
      REAL*8 E_nu
      REAL*8 xsect_total

C  Local Variable Declarations

      REAL*8 energy, xsect(energy_size)
      INTEGER e_dum, e_index

      INTEGER xs_data_entries
      PARAMETER (xs_data_entries = xs_nu_208Pb_cc_Engel_2n_entries)  ! from 'xsect_size'
      CHARACTER*40 xs_data_file

      REAL*8 E_thres

      REAL*8 E_data(xs_data_entries), 
     +       xs_nu_data(xs_data_entries), xs_nubar_data(xs_data_entries)
      REAL*8 b_nu_xs(aspline_size), c_nu_xs(aspline_size),
     +       d_nu_xs(aspline_size)
      REAL*8 b_nubar_xs(aspline_size), c_nubar_xs(aspline_size),
     +       d_nubar_xs(aspline_size)
      REAL*8 xsect_fit

      INTEGER n_calls, last_interaction

      INTEGER i

      SAVE n_calls, last_interaction

      SAVE E_data, xs_nu_data, xs_nubar_data
      SAVE b_nu_xs, c_nu_xs, d_nu_xs
      SAVE b_nubar_xs, c_nubar_xs, d_nubar_xs
      SAVE xsect

************************************************************************
*  Initializations
************************************************************************

      DATA n_calls             /0/
      DATA last_interaction    /0/

      DATA E_data              /xs_data_entries * 0./,
     +     xs_nu_data          /xs_data_entries * 0./,
     +     xs_nubar_data       /xs_data_entries * 0./

      DATA b_nu_xs             /aspline_size * 0./,
     +     c_nu_xs             /aspline_size * 0./,
     +     d_nu_xs             /aspline_size * 0./,
     +     b_nubar_xs          /aspline_size * 0./,
     +     c_nubar_xs          /aspline_size * 0./,
     +     d_nubar_xs          /aspline_size * 0./

      DATA xsect               /energy_size * 0./

      xs_data_file = 'dat/xsect/xs_nu_208Pb_cc_Engel_2n.dat'

*************************************************************************

* reset counter for the next interaction

      IF (interaction .NE. last_interaction) n_calls = 0

* increment counter for this interaction

      n_calls = n_calls + 1
      
* if we've gone through the energy loop once, just use saved values
* (doesn't work if we bypass the integrations over time and energy by
*  reading in the count sum totals, therefore always calculate the cross
*  section explicitly in that case)
      
      IF (      n_calls .GT. energy_entries
     +    .AND. gen_flag .NE. 3) THEN

         xsect_total = xsect(this_energy_index)

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

         CALL get_threshold(interaction,E_thres)

**************************************************************************
*  Get cross section data to spline fit, convert data to logarithms because
*  they are easier to fit, then evaluate the Akima spline fit coefficients 
*  (done once per interaction)
**************************************************************************

         IF (n_calls .EQ. 1) THEN

c            WRITE (*,*) '  --> xs_nu_208Pb_nc.f: filling data arrays...' 

            OPEN (unit = 50, file = xs_data_file, status = 'old')

            DO 20 i = 1, xs_data_entries
               READ (50,*) E_data(i), xs_nu_data(i)

               xs_nu_data(i) = dlog10(xs_nu_data(i))

20          CONTINUE

            CLOSE (unit = 50)

            IF (interaction .EQ. 80) THEN

               CALL aspline(xs_data_entries,E_data,xs_nu_data,
     +                      b_nu_xs,c_nu_xs,d_nu_xs)

C            ELSEIF (     interaction .EQ. 85) THEN

C               CALL aspline(xs_data_entries,E_data,xs_nubar_data,
C     +                      b_nubar_xs,c_nubar_xs,d_nubar_xs)

            ENDIF
         ENDIF     ! condition for first call

*****************************************************************************
*  Evaluate the Akima spline fit (done for every neutrino for each
*  interaction)
*****************************************************************************

         IF (     interaction .EQ. 80) THEN

            CALL aspline_eval(xs_data_entries,E_data,xs_nu_data,
     +                        b_nu_xs,c_nu_xs,d_nu_xs,E_nu,xsect_fit)

C         ELSEIF (        interaction .EQ. 85) THEN

C            CALL aspline_eval(xs_data_entries,E_data,xs_nubar_data,
C     +                        b_nubar_xs,c_nubar_xs,d_nubar_xs,
C     +                        E_nu,xsect_fit)

         ENDIF

* cross section in [m^2]

         IF (E_nu .GT. E_thres) THEN
            xsect_total = 10.d0**(xsect_fit)
            xsect_total = xsect_total * cm2m**2
         ELSE
            xsect_total = 0.d0
         ENDIF

* fill cross section array since we haven't gone through the entire
* energy spectrum yet
 
C         WRITE (*,*) E_nu, xsect_total
           
         xsect(this_energy_index) = xsect_total

      ENDIF   ! calculate or use saved array values

* remember some things about this time through the routine
 
      last_interaction = interaction

      RETURN
      END
