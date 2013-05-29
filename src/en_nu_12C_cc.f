**************************************************************************
*   Filename: en_nu_12C_cc.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_12C_cc(interaction,E_nu,E_part,decay_delay)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION E_thres
      REAL*8 E_part, E_nu, decay_delay
      INTEGER n_particle, interaction

      DATA n_particle /1/

      SAVE n_particle

*****************************************************************************

      decay_delay = 0
      CALL get_threshold(interaction,E_thres)

      IF (n_particle .EQ. 1) THEN

          E_part = E_nu - E_thres 

          n_particle = 2

      ELSEIF (n_particle .EQ. 2) THEN

          E_part = E_thres

          n_particle = 1
        
          IF (     interaction .EQ. 59
     +        .OR. interaction .EQ.  4) THEN
                decay_delay = 1100
          ELSEIF (     interaction .EQ. 64
     +             .OR. interaction .EQ. 16) THEN
                decay_delay = 2020
          ENDIF

      ENDIF

      RETURN
      END
