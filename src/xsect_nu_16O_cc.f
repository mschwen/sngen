***********************************************************************
*  Filename: xsect_nu_16O_cc.f
***********************************************************************
*
*  This routine returns the total nu + 16O CC cross section for a given
*  neutrino energy.  The cross sections used are those taken from:
*       o Kuramoto et al.  (1989)
*
* Also included are the cross sections that Burrows used:
*       o Burrows et al.  (1992)
*  This routine uses the parameterization provided by Adam Burrows.
*
*  Inputs:
*       o  interaction id number (not array)
*       o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy (m^2)
*
***********************************************************************
     
      SUBROUTINE xsect_nu_16O_cc(interaction,E_nu,xsect_total)

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

C  Global Variable Declarations

      INTEGER          interaction
      REAL*8           E_nu
      REAL*8           xsect_total

C  Local Variable Declarations

      LOGICAL          first_call
      REAL*8           energy, xsect(energy_size)
      INTEGER          e_dum, e_index

      INTEGER          xs_data_entries
      PARAMETER        (xs_data_entries = xs_nu_O_entries)  ! default
      CHARACTER*40     xs_data_file

      REAL*8           E_thres

      REAL*8           E_data(xs_data_entries), 
     +                 xs_nu_data(xs_data_entries), 
     +                 xs_nubar_data(xs_data_entries)
      REAL*8           b_nu_xs(aspline_size), 
     +                 c_nu_xs(aspline_size),
     +                 d_nu_xs(aspline_size)
      REAL*8           b_nubar_xs(aspline_size), 
     +                 c_nubar_xs(aspline_size),
     +                 d_nubar_xs(aspline_size)
      REAL*8           xsect_fit

      INTEGER          n_calls, last_interaction

      INTEGER          i
 
      SAVE             first_call
      SAVE             n_calls, last_interaction
      SAVE             E_data, xs_nu_data, xs_nubar_data,
     +                 b_nu_xs, c_nu_xs, d_nu_xs,
     +                 b_nubar_xs, c_nubar_xs, d_nubar_xs
      SAVE             xsect

************************************************************************
*  Initializations
************************************************************************

      DATA             first_call       /.true./
      DATA             n_calls          /0/
      DATA             last_interaction /0/
      
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

      xs_data_file = 'dat/xsect/xs_nu_O_cc.dat'           ! default

*************************************************************************

* reset counter for the next interaction
      
      IF (interaction .NE. last_interaction) n_calls = 0

* increment counter for this interaction

      n_calls = n_calls + 1

* if we've gone through the energy loop once, just use saved values
* (doesn't work if we bypass the integrations over time and energy by
*  reading in the count sum totals, therefore always calculate the cross
*  section explicitly in that case)
      
      IF (      (n_calls .GT. energy_entries)
     +    .AND. (gen_flag .NE. 3)          ) THEN

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

         IF (E_nu .GT. E_thres) THEN

            IF (      sn_model_switch .EQ. 1
     +          .AND. use_author_params .EQ. 1) THEN

*************************************************************************
*  Burrows cross section (courtesy Adam Burrows)
*  Note: Burrows imposes a 5 MeV threshold on all SNO reactions,
*        but in this case, nu_e + 16O CC efficiency only done with energy 
*        threshold for nu + 18O reaction
*************************************************************************

               IF (     interaction .EQ. 5
     +             .OR. interaction .EQ. 6
     +             .OR. interaction .EQ. 7) THEN

                  IF (E_nu .LT. 10.d0) THEN
                     xsect_total = 0.0259d0 * (E_nu - 1.66d0)**2.39d0
                  ELSEIF (E_nu .LT. 22.05d0) THEN
                     xsect_total = 0.0764d0 * (E_nu - 1.66d0)**1.88d0
                  ELSEIF (E_nu .LT. 50.0d0) THEN
                     xsect_total = 0.0556d0 * (E_nu - 15.4d0)**3.16d0
                  ELSE
                     xsect_total = 0.00596d0 * (E_nu - 15.4d0)**3.79d0
                  ENDIF

               ELSEIF (     interaction .EQ. 17
     +                 .OR. interaction .EQ. 18) THEN
   
                  IF (E_nu .LT. 20.0d0) THEN
                     xsect_total = 0.0212d0 * (E_nu - 11.4d0)**3.01d0
                  ELSEIF (E_nu .LT. 29.5d0) THEN
                     xsect_total = 0.00452d0 * (E_nu - 11.4d0)**3.73d0
                  ELSEIF (E_nu .LT. 41.1d0) THEN
                     xsect_total = 0.0121d0 * (E_nu - 11.4d0)**3.39d0
                  ELSE
                     xsect_total = 0.0410d0 * (E_nu - 11.4d0)**3.03d0
                  ENDIF
               ENDIF  ! check for interaction type
* cross section in [m^2]

               xsect_total = xsect_total * 1.d-44 * cm2m**2

            ELSE    ! use the default method

**************************************************************************
*  Get cross section data to spline fit, convert data to logarithms because
*  they are easier to fit, then evaluate the Akima spline fit coefficients
*  (done once so do both nu_e and nu_e_bar at the same time)
**************************************************************************

               IF (first_call) THEN

                  first_call = .false.
 
                  OPEN (unit = 50, file = xs_data_file, status = 'old')
      
                  DO 20 i = 1, xs_data_entries
                     READ (50,*) E_data(i), 
     +                           xs_nu_data(i), xs_nubar_data(i)
            
                     xs_nu_data(i)    = dlog10(xs_nu_data(i))
                     xs_nubar_data(i) = dlog10(xs_nubar_data(i))

20                CONTINUE

                  CLOSE (unit = 50)

* interaction 5,6,7
                  CALL aspline(xs_data_entries,E_data,xs_nu_data,
     +                         b_nu_xs,c_nu_xs,d_nu_xs)
               
* interaction 17,18
                  CALL aspline(xs_data_entries,E_data,xs_nubar_data,
     +                         b_nubar_xs,c_nubar_xs,d_nubar_xs)

               ENDIF     ! condition for first call
     
*****************************************************************************
*  Evaluate the Akima spline fit (done for every neutrino for each
*  interaction)
*****************************************************************************
               
               IF (     interaction .EQ. 5
     +             .OR. interaction .EQ. 6
     +             .OR. interaction .EQ. 7) THEN

                  CALL aspline_eval(xs_data_entries,
     +                              E_data,xs_nu_data,
     +                              b_nu_xs,c_nu_xs,d_nu_xs,
     +                              E_nu,xsect_fit)

               ELSEIF (     interaction .EQ. 17
     +                 .OR. interaction .EQ. 18) THEN

                  CALL aspline_eval(xs_data_entries,
     +                              E_data,xs_nubar_data,
     +                              b_nubar_xs,c_nubar_xs,d_nubar_xs,
     +                              E_nu,xsect_fit)

               ENDIF  ! check for which interaction

* cross section in [m^2]
     
               xsect_total = 10.d0**(xsect_fit)
               xsect_total = xsect_total * cm2m**2

               IF (      sn_model_switch .EQ. 2         ! Beacom & Vogel
     +             .AND. use_author_params .EQ. 1) THEN

                  IF ((E_nu - E_thres) .LT. 5.d0) xsect_total = 0.d0 ! det thres

               ENDIF ! check for use of BeacomVogel author parameters
            ENDIF    ! check for use of Burrows author parameters
         ELSE
            xsect_total = 0.d0
         ENDIF     ! check for E_nu > E_thres
     
* fill cross section array since we haven't gone through the entire
* energy spectrum yet

         xsect(this_energy_index) = xsect_total

c         WRITE (*,*) 'E_nu, xsect_total =', E_nu, xsect_total
         
      ENDIF   ! calculate or use saved array values

* remember some things about this time through the routine
         
      last_interaction = interaction
            
      RETURN 
      END
