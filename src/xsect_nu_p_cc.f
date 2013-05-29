***********************************************************************
*  Filename: xsect_nu_p_cc.f
***********************************************************************
*
*  This routine returns the total nu_e_bar + p CC cross section for a
*  given neutrino energy.  Two sets of cross sections are included:
*	o Llewellyn Smith with modifications suggested by Vogel & Beacom
*	o Burrows et al. (1992)
*  Comment out the appropriate lines in the code below to change the
*  default.
*
*  For more detail on the Llewellyn Smith cross section, see the routine
*  'xsect_nu_p_diff.f' where the formalism is defined.
*
*  Inputs:
*	o  interaction id number (not array)
*	o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy (m^2)
*
***********************************************************************
     
      SUBROUTINE xsect_nu_p_cc(interaction,E_nu,xsect_total)

      IMPLICIT none

      INCLUDE 'energy_size.inc'
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'
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

      REAL*8           energy, xsect(energy_size)
      INTEGER          e_dum, e_index

* variables for Vogel & Beacom calculations
      INTEGER          i, istep
      PARAMETER        (istep = 600)          ! number of integration steps
      DOUBLE PRECISION a, b, c,
     +                 E_e_min, E_e_max,
     +                 q2, q2_min, q2_max, q2_step,
     +                 sum_xsect,
     +                 xsect_nu_p_diff         ! function

      DOUBLE PRECISION m_e, M, m_n, m_p
      PARAMETER        (m_e = m_electron)
      PARAMETER        (m_n = m_neutron)
      PARAMETER        (m_p = m_proton)
      PARAMETER        (M = m_nucleon)

      REAL*8           E_thres

      INTEGER          n_calls, last_interaction

      SAVE             n_calls, last_interaction,
     +                 xsect

********************************************************************
*  Initializations
********************************************************************
      
      DATA n_calls           /0/
      DATA last_interaction  /0/
      DATA xsect             /energy_size * 0./
      
************************************************************************

* reset counter for the next interaction
     
      IF (interaction .NE. last_interaction) n_calls = 0

* increment counter for this interaction

      n_calls = n_calls + 1
      
* if we've gone through the energy loop once, just use saved values
* (doesn't work if we bypass the integrations over time and energy by
*  reading in the count sum totals, therefore always calculate the cross
*  section explicitly in that case)
      
      IF (      (n_calls .GT. energy_entries)
     +    .AND. (gen_flag .NE. 3)            ) THEN

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

* compute cross section  [m^2]

         IF (E_nu .GT. E_thres) THEN
            IF (      sn_model_switch .EQ. 1
     +          .AND. use_author_params .EQ. 1) THEN

*************************************************************************
*  Burrows cross section (courtesy Adam Burrows)
*  (doesn't reproduce paper results, though;  may have been updated)
*  Note: Burrows imposes a 5 MeV threshold on all SNO reactions
*************************************************************************

               IF ((E_nu - E_thres) .GE. 5.d0) THEN
                  xsect_total = 9.41d0 * (E_nu)**2 
     +                         * (1.d0 - 1.29335d0 / E_nu)
     +                         * dsqrt((1.d0 - 1.29335d0/E_nu)**2
     +                          - (0.511d0/E_nu)**2)
     +                         * (1.d0 - 0.00325d0*(E_nu - 0.64667d0))

* Burrows' generator gets 325.450439 counts with the cross section above,
* but his paper claims 331, so scale to get agreement

                  xsect_total = 1.01705d0 * xsect_total ! reproduce paper result

                  xsect_total = xsect_total * 1.d-44 * cm2m**2
               ELSE
                  xsect_total = 0.d0
               ENDIF  ! not above detector threshold Burrows imposes

            ELSE      ! use the default method with Burrows' model

**************************************************************************
*  Llewellyn-Smith (with modifications from Vogel & Beacom) method:
*  Integrate the differential quasi-elastic cross section.
**************************************************************************

* determine the range of E_e allowed from the quadratic equation

         a = (2.d0 * (E_nu + m_p))**2 - (2.d0 * E_nu)**2
         b = 2.d0 * (2.d0 * M * delta_m - 2.d0 * m_p * E_nu - m_e**2) 
     +      * (2.d0 * (E_nu + m_p))
         c = (2.d0 * M * delta_m - 2.d0 * m_p * E_nu - m_e**2)**2 
     +      + (2.d0 * E_nu * m_e)**2

         E_e_min = (-b - dsqrt(b**2 - 4.d0*a*c))/(2.d0*a)
         E_e_max = (-b + dsqrt(b**2 - 4.d0*a*c))/(2.d0*a)

* convert to q^2 range (Vogel & Beacom Equation 6)

         q2_min = 2.d0*M*delta_m - 2.d0*m_p*(E_nu - E_e_min)
         q2_max = 2.d0*M*delta_m - 2.d0*m_p*(E_nu - E_e_max)

* setup integration steps and integrate over q^2

         q2_step = (q2_max - q2_min) / dble(istep)

         sum_xsect = 0.d0

         DO i = 1, istep
            q2 = q2_min + (dble(i) - 0.5d0)*q2_step
            sum_xsect = sum_xsect + xsect_nu_p_diff(E_nu,q2) * q2_step
         END DO

         xsect_total = sum_xsect * cm2m**2

               IF (      sn_model_switch .EQ. 2         ! Beacom & Vogel
     +             .AND. use_author_params .EQ. 1) THEN
                  xsect_total = 0.991913746d0 * xsect_total  ! old g_A value
                  xsect_total = 1.0682d0 * xsect_total  ! convert to Vogel
                 IF ((E_nu - E_thres) .LT. 5.d0) xsect_total = 0.d0 ! det thres
               ENDIF  ! use Beacom & Vogel parameters
            ENDIF     ! check for using Burrows author parameters
         ELSE
            xsect_total = 0.d0
         ENDIF

* fill cross section array since we haven't gone through the entire 
* energy spectrum yet

         xsect(this_energy_index) = xsect_total

      ENDIF   ! calculate or use saved array values

* remember some things about this time through the routine
      
      last_interaction = interaction

      RETURN
      END
