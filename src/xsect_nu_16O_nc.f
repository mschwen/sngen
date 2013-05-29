***********************************************************************
*  Filename: xsect_nu_16O_nc.f
***********************************************************************
*
*  This routine returns the total NC cross section for the reactions
*  involving 'nu + 160' NC interactions for a given neutrino energy.
*
*  These cross section ratios are adapted from the work of:
*        Langanke, Vogel & Kolbe PRL 76 2629 (1996)
*  
*  The parameterization for the total cross section is taken from:
*        Beacom & Vogel hep-ph/9802424 v2 10 June 1998
*     
*  Inputs:
*       o  interaction
*	o  neutrino energy (MeV)
*  Outputs:
*       o  total interaction cross section for the given energy
*
***********************************************************************
     
      SUBROUTINE xsect_nu_16O_nc(interaction,E_nu,xsect_total)

      IMPLICIT none

      INCLUDE 'energy_size.inc'
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

      REAL*8 E_thres

      REAL*8 sig_0, branch_frac

* [Beacom & Vogel hep-ph/9802424 v2 10 June 1998]
* (made more general by removing gamma final state branching fraction)
* -> 0.75e-47 given in hep-ph/9802424 v2
* -> other branching fractions given in PRL 76 2629 (1996) 
*
*  sig_0 = 0.75e-47 / branch1 / branch2
*
*   --> branch1 = sum(p or n) / sum(total)
*               = (3.75 + 1.76 + 2.93 + 1.29) / (5.90 + 4.48) 
*               = 0.937
*
*   --> branch2 = sum(n or p with gamma) / sum(p or n)
*               = (1.41 + 0.37 + 1.09 + 0.28) / (3.75 + 1.76 + 2.93 + 1.29)
*               = 0.324
*
*     ==> sig_0 = 0.75e-47 / 0.937 / 0.324
*               = 2.47e-47

      PARAMETER (sig_0 = 2.471428571d-47)

      INTEGER n_calls, last_interaction

      SAVE n_calls, last_interaction

************************************************************************
*  Initializations
************************************************************************

      DATA n_calls             /0/
      DATA last_interaction    /0/

      DATA xsect               /energy_size * 0./

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

         CALL get_threshold (interaction,E_thres)

         IF (E_nu .GT. E_thres) THEN

**************************************************************************
*  Assign branching fractions for the various reactions
**************************************************************************

*** nu_mu, nu_tau

* n
            IF (     interaction .EQ. 25
     +          .OR. interaction .EQ. 43) THEN

               branch_frac = 1.76d0 / 5.90d0

* p + gamma
            ELSEIF (     interaction .EQ. 26
     +              .OR. interaction .EQ. 27
     +              .OR. interaction .EQ. 44
     +              .OR. interaction .EQ. 45) THEN

               branch_frac = 1.41d0 / 5.90d0

* n + gamma
            ELSEIF (     interaction .EQ. 28
     +              .OR. interaction .EQ. 29
     +              .OR. interaction .EQ. 46
     +              .OR. interaction .EQ. 47) THEN

               branch_frac = 0.37d0 / 5.90d0

*** nu_mubar, nu_taubar

* n
            ELSEIF (     interaction .EQ. 34
     +              .OR. interaction .EQ. 52) THEN

               branch_frac = 1.29d0 / 4.48d0

* p + gamma
            ELSEIF (     interaction .EQ. 35
     +              .OR. interaction .EQ. 36
     +              .OR. interaction .EQ. 53
     +              .OR. interaction .EQ. 54) THEN

               branch_frac = 1.09d0 / 4.48d0

* n + gamma
            ELSEIF (     interaction .EQ. 37
     +              .OR. interaction .EQ. 38
     +              .OR. interaction .EQ. 55
     +              .OR. interaction .EQ. 56) THEN

               branch_frac = 0.28d0 / 4.48d0

            ENDIF

            branch_frac = 0.5d0 * branch_frac  ! distinguish nu_mu, nu_tau

            xsect_total = sig_0 * branch_frac * (E_nu - E_thres)**4

* convert cross section from [cm^2] to [m^2]

            xsect_total = xsect_total * cm2m**2

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
