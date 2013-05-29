***********************************************************************
*  Filename: xsect_nu_p_es.f
***********************************************************************
*
*  This routine returns the total nu + p ES cross section for a given
*  neutrino energy.
*     
*  Inputs:
*       o  interaction id  (not array)
*	o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy (m^2)
*
***********************************************************************
     
      SUBROUTINE xsect_nu_p_es(interaction,E_nu,xsect_total)

      IMPLICIT none

      INCLUDE 'energy_size.inc'
      INCLUDE 'id_nu.inc'
      INCLUDE 'physics_constants.inc'
      INCLUDE 'trig_constants.inc'
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
 
      COMMON           /det_threshold/
     +                 E_min
      REAL*8           E_min

C  Global Variable Declarations

      INTEGER          interaction  ! input
      REAL*8           E_nu         ! input
      REAL*8           xsect_total  ! output   

C  Local Variable Declarations

      REAL*8           energy, xsect(energy_size) 
      INTEGER          e_dum, e_index

      INTEGER          nu_type
      REAL             threshold

C  Cross section Calculation

      LOGICAL          first_call

      DOUBLE PRECISION T_p_min, T_p_max
      DOUBLE PRECISION term1, term2, term3
      DOUBLE PRECISION sig_p
      DOUBLE PRECISION cv, ca

      INTEGER          n_calls, last_interaction
      
      SAVE             n_calls, last_interaction
      SAVE             xsect

*********************************************************************
*  Initializations
*********************************************************************

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

c         write (*,*) '  --> using saved xsect values!'
c         write (*,*) '      n_calls, energy_entries=',
c     +                n_calls, energy_entries

         xsect_total = xsect(this_energy_index)

c         write (*,*) 'this_energy_index, xsect_total=',
c     +                this_energy_index, xsect_total

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

         CALL get_nu_type(interaction,nu_type)

********************************************************************
* Cross section calculation
********************************************************************

            IF (     nu_type .EQ. id_nu_e
     +          .OR. nu_type .EQ. id_nu_mu    
     +          .OR. nu_type .EQ. id_nu_tau) THEN

                cv = 0.04
                ca = 1.27 / 2

            ELSE

                cv = 0.04
                ca = -1.27 / 2

            ENDIF


**************************************************************************
*  Cross section energy dependent stuff (including T_min dependence)
**************************************************************************

* maximum electron kinetic energy
               
               T_p_max = 2.d0 * E_nu**2 / (2.d0 * E_nu + m_proton)
               T_p_min = E_min               

* various cross section terms

               term1 = (cv + ca)**2 * (T_p_max - T_p_min) * E_nu**2
               term2 = (cv - ca)**2 * (E_nu**2 * (T_p_max - T_p_min)
     +                 - E_nu*(T_p_max**2 - T_p_min**2) + (T_p_max**3
     +                 - T_p_min**3) / 3.d0)
               term3 = - (cv**2 - ca**2) * m_proton * (T_p_max**2 - 
     +                 T_p_min**2) / 2.d0

               sig_p = (G_fermi**2 * m_proton * hbar_c2) / 
     +                 (2*pi*(E_nu)**2)

* add all the terms and factor in the overall cross section constant

               xsect_total = sig_p * (term1 + term2 + term3)

C               sig_p = sig_p * E_nu**2

C                IF (E_nu .EQ. 170) THEN
C                WRITE (*,*) term1, term2, term3, T_p_max, sig_p
C                ENDIF

* convert cross section from [cm] to [m]

            xsect_total = xsect_total * cm2m**2

* using detector threshold can give negative values, so set to 0 in that case

            IF (xsect_total .LT. 0.d0) xsect_total = 0.d0

            CALL get_threshold(interaction, threshold)

            IF (E_nu .LT. threshold) xsect_total = 0.d0 

* fill cross section array since we haven't gone through the entire 
* energy spectrum yet

         xsect(this_energy_index) = xsect_total

c         WRITE (*,*) 'this_energy_index, E_nu, xsect_total =',
c     +                this_energy_index, E_nu, xsect_total

      ENDIF  ! use saved values or do calculation

* remember some things about this time through the routine

      last_interaction = interaction

      RETURN
      END
