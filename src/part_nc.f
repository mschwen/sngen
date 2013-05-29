***************************************************************************
* Filename: part_nc.f
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

      SUBROUTINE part_nc(interaction,E_nu,id_part,E_part,
     +                   sn_dir_cos,part_dir)
   
      IMPLICIT none

      INCLUDE 'id_part.inc'

C   Global Variable Declarations

      INTEGER interaction, id_part
      REAL*8 E_nu, E_part
      REAL*8 sn_dir_cos(*), part_dir(3)

***************************************************************************

*** nu + d  NC reactions

      IF (     interaction .EQ. 8
     +    .OR. interaction .EQ. 19
     +    .OR. interaction .EQ. 23
     +    .OR. interaction .EQ. 32
     +    .OR. interaction .EQ. 41
     +    .OR. interaction .EQ. 50) THEN

         IF (id_part .EQ. id_neutron) THEN
            CALL en_nu_d_nc(interaction,E_nu,E_part)
            CALL dir_nu_d_nc(sn_dir_cos,part_dir)
         ENDIF

*** HALO 2 neutron NC reactions

      ELSEIF (     interaction .EQ. 82
     +        .OR. interaction .EQ. 87
     +        .OR. interaction .EQ. 90
     +        .OR. interaction .EQ. 93
     +        .OR. interaction .EQ. 96
     +        .OR. interaction .EQ. 99) THEN

         CALL dir_ran(sn_dir_cos,part_dir)
         CALL en_nu_Pb_nc_2n(interaction,E_nu,E_part)

*** HALO 1 neutron NC reactions

      ELSEIF (     interaction .EQ. 81
     +        .OR. interaction .EQ. 86
     +        .OR. interaction .EQ. 89
     +        .OR. interaction .EQ. 92
     +        .OR. interaction .EQ. 95
     +        .OR. interaction .EQ. 98) THEN

         CALL dir_ran(sn_dir_cos,part_dir)
         CALL en_nu_Pb_nc_1n(interaction,E_nu,E_part)

*** LAB NC reactions

      ELSEIF (     interaction .EQ. 60
     +        .OR. interaction .EQ. 65
     +        .OR. interaction .EQ. 68
     +        .OR. interaction .EQ. 71
     +        .OR. interaction .EQ. 74
     +        .OR. interaction .EQ. 77) THEN

         CALL dir_ran(sn_dir_cos,part_dir)
         CALL en_nu_12C_nc(E_part)

*** AV 12C NC reactions

      ELSEIF (     interaction .EQ. 9
     +        .OR. interaction .EQ. 20
     +        .OR. interaction .EQ. 24
     +        .OR. interaction .EQ. 33
     +        .OR. interaction .EQ. 42
     +        .OR. interaction .EQ. 51) THEN

         CALL en_nu_12C_nc(E_part)
         CALL dir_ran(sn_dir_cos,part_dir)         

*** H20 O16 NC -> 15O reactions

       ELSEIF (     interaction .EQ. 28
     +         .OR. interaction .EQ. 37
     +         .OR. interaction .EQ. 46
     +         .OR. interaction .EQ. 55) THEN

          CALL dir_ran(sn_dir_cos,part_dir)
          CALL en_nu_16O_nc_15O(E_part)

*** H20 O16 NC -> 15N reactions

       ELSEIF (     interaction .EQ. 26
     +         .OR. interaction .EQ. 35
     +         .OR. interaction .EQ. 44
     +         .OR. interaction .EQ. 53) THEN

          CALL dir_ran(sn_dir_cos,part_dir)
          CALL en_nu_16O_nc_15N(E_part)

       ENDIF

      RETURN
      END
