***************************************************************************
* Filename: part_cc.f
***************************************************************************
*
*   This routine simply redirects the interaction information to the
*   appropriate reaction routine 
*
*   Inputs:
*	o  interaction id number (not array)
*	o  delayed neutrino interaction time (ms)
*	o  neutrino energy (MeV)
*   Outputs:
*	o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE part_cc(interaction,E_nu,id_part,E_part,
     +                   sn_dir_cos,part_dir,redo,decay_delay)
   
      IMPLICIT none

      INCLUDE 'id_part.inc'

C   Global Variable Declarations

      INTEGER interaction, id_part
      REAL*8 E_nu, E_part, decay_delay
      REAL*8 sn_dir_cos(*), part_dir(3)
      INTEGER redo

***************************************************************************

      redo = 0

***************************************************************************

*** nu_e + d --> p + p + e^- CC reaction

      IF (interaction .EQ. 3) THEN
         CALL en_nu_d_cc(interaction,E_nu,E_part)
         CALL dir_nu_d_cc(E_part,sn_dir_cos,part_dir)

*** nu_e_bar + d --> n + n + e^+ CC reaction

      ELSEIF (interaction .EQ. 15) THEN
         IF (id_part .EQ. id_positron) THEN
            CALL en_nu_d_cc(interaction,E_nu,E_part)
            CALL dir_nu_d_cc(E_part,sn_dir_cos,part_dir)
         ELSEIF (id_part .EQ. id_neutron) THEN
            CALL en_nu_d_n_cc(E_nu,E_part,redo)
            IF (redo .EQ. 1) RETURN
            CALL dir_nu_d_n_cc(E_nu,part_dir)
         ENDIF

*** nu_e_bar + p --> n + e^+ CC reactions

      ELSEIF (     interaction .EQ. 12
     +        .OR. interaction .EQ. 13
     +        .OR. interaction .EQ. 14
     +        .OR. interaction .EQ. 63) THEN
         CALL en_nu_p_cc(E_nu,id_part,E_part)
         CALL dir_nu_p_cc(E_nu,id_part,E_part,sn_dir_cos,part_dir)

*** AV 12C CC reactions

      ELSEIF (     interaction .EQ. 4
     +        .OR. interaction .EQ. 16) THEN
         CALL en_nu_12C_cc(interaction,E_nu,E_part,decay_delay)
         CALL dir_ran(sn_dir_cos,part_dir)

*** LAB 12C CC reactions

      ELSEIF (     interaction .EQ. 59
     +        .OR. interaction .EQ. 64) THEN
         CALL en_nu_12C_cc(interaction,E_nu,E_part,decay_delay)
         CALL dir_ran(sn_dir_cos,part_dir)

*** H20 16O CC reactions

      ELSEIF (     interaction .EQ. 5
     +        .OR. interaction .EQ. 17) THEN
         CALL en_nu_16O_cc(interaction,E_nu,E_part)
         CALL dir_ran(sn_dir_cos,part_dir)

*** HALO 1 neutron CC reactions

      ELSEIF (     interaction .EQ. 79
     +        .OR. interaction .EQ. 84) THEN
          CALL en_nu_Pb_cc_1n(interaction,E_nu,E_part)
          CALL dir_ran(sn_dir_cos,part_dir)

*** HALO 2 neutron CC reactions

      ELSEIF (     interaction .EQ. 80
     +        .OR. interaction .EQ. 85) THEN
          CALL en_nu_Pb_cc_2n(interaction,E_nu,E_part)
          CALL dir_ran(sn_dir_cos,part_dir)

      ENDIF

      RETURN
      END
