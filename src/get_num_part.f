***********************************************************************
*  Filename: get_num_part.f
***********************************************************************
*
*  This routine returns the number of final state detected particles
*  as well as their names (id numbers).  See also 'id_part.inc'.
*   
*  Inputs:
*       o  interaction id number (not array)
*  Outputs:
*       o  number of final state detected particles (integer)
*	o  particle names (array of ids)
*
***********************************************************************
     
      SUBROUTINE get_num_part(interaction,npart,id_part)

      IMPLICIT none

      INCLUDE 'id_int.inc'
      INCLUDE 'id_part.inc'

C  Global Variable Declarations

      INTEGER interaction
      INTEGER npart
      INTEGER id_part(3)      ! maximum 3 final state detected particles

C  Local Variable Declarations

      INTEGER int_type
      INTEGER i
      REAL random

************************************************************************
*  Initializations
************************************************************************

      npart = 0

      DO 10 i = 1, 3
         id_part(i) = 0
10    CONTINUE

************************************************************************

      CALL get_int_type(interaction,int_type)

*** elastic scattering reactions

      IF (int_type .EQ. id_nu_es) THEN
         IF (     interaction .EQ. 58
     +       .OR. interaction .EQ. 62
     +       .OR. interaction .EQ. 67
     +       .OR. interaction .EQ. 70
     +       .OR. interaction .EQ. 73
     +       .OR. interaction .EQ. 76) THEN

             npart = 1
             id_part(1) = id_proton
         ELSE 
             npart = 1
             id_part(1) = id_electron
         ENDIF

*** neutral current reactions

      ELSEIF (int_type .EQ. id_nu_nc) THEN
         IF (     interaction .EQ. 8
     +       .OR. interaction .EQ. 19
     +       .OR. interaction .EQ. 23
     +       .OR. interaction .EQ. 25
     +       .OR. interaction .EQ. 32
     +       .OR. interaction .EQ. 34
     +       .OR. interaction .EQ. 41
     +       .OR. interaction .EQ. 43
     +       .OR. interaction .EQ. 50
     +       .OR. interaction .EQ. 52) THEN
            npart = 1
            id_part(1) = id_neutron

         ELSEIF (     interaction .EQ. 9
     +           .OR. interaction .EQ. 20
     +           .OR. interaction .EQ. 24
     +           .OR. interaction .EQ. 33
     +           .OR. interaction .EQ. 42
     +           .OR. interaction .EQ. 51
     +           .OR. interaction .EQ. 60
     +           .OR. interaction .EQ. 65
     +           .OR. interaction .EQ. 68
     +           .OR. interaction .EQ. 71
     +           .OR. interaction .EQ. 74
     +           .OR. interaction .EQ. 77) THEN
                npart = 3
                id_part(1) = id_gamma
                id_part(2) = id_gamma
                id_part(3) = id_gamma
             
         ELSEIF (     interaction .EQ. 28
     +           .OR. interaction .EQ. 29
     +           .OR. interaction .EQ. 38
     +           .OR. interaction .EQ. 45
     +           .OR. interaction .EQ. 56) THEN
            npart = 2
            id_part(1) = id_neutron	
            id_part(2) = id_gamma

         ELSEIF (     interaction .EQ. 81
     +           .OR. interaction .EQ. 86
     +           .OR. interaction .EQ. 89
     +           .OR. interaction .EQ. 92
     +           .OR. interaction .EQ. 95
     +           .OR. interaction .EQ. 98) THEN
            npart = 1
            id_part(1) = id_neutron

         ELSEIF (     interaction .EQ. 82
     +           .OR. interaction .EQ. 87
     +           .OR. interaction .EQ. 90
     +           .OR. interaction .EQ. 93
     +           .OR. interaction .EQ. 96
     +           .OR. interaction .EQ. 99) THEN
            npart = 2
            id_part(1) = id_neutron
            id_part(2) = id_neutron

         ELSEIF (     interaction .EQ. 35
     +           .OR. interaction .EQ. 44
     +           .OR. interaction .EQ. 53) THEN
            npart = 1
            id_part(1) = id_gamma

         ELSEIF (     interaction .EQ. 37
     +           .OR. interaction .EQ. 46
     +           .OR. interaction .EQ. 55) THEN
            npart = 2
            id_part(1) = id_gamma
            id_part(2) = id_gamma
         ENDIF

         
***  charged current reactions

      ELSE
         IF (     interaction .EQ. 3
     +       .OR. interaction .EQ. 5
     +       .OR. interaction .EQ. 6
     +       .OR. interaction .EQ. 7) THEN
            npart = 1
            id_part(1) = id_electron

         ELSEIF (interaction .EQ. 12) THEN
            npart = 1
            id_part(1) = id_positron

         ELSEIF (     interaction .EQ. 13
     +           .OR. interaction .EQ. 14
     +           .OR. interaction .EQ. 63) THEN
            npart = 2
            id_part(1) = id_positron
            id_part(2) = id_neutron

         ELSEIF (interaction .EQ. 15) THEN
            npart = 3
            id_part(1) = id_positron
            id_part(2) = id_neutron
            id_part(3) = id_neutron

         ELSEIF (     interaction .EQ. 16
     +           .OR. interaction .EQ. 17
     +           .OR. interaction .EQ. 18) THEN
            npart = 1
            id_part(1) = id_positron

         ELSEIF (interaction .EQ. 79) THEN
            npart = 2
            id_part(1) = id_neutron
            id_part(2) = id_electron

         ELSEIF (interaction .EQ. 80) THEN
            npart = 3
            id_part(1) = id_neutron
            id_part(2) = id_neutron
            id_part(3) = id_electron 

         ELSEIF (     interaction .EQ. 4
     +           .OR. interaction .EQ. 59) THEN
            npart = 2
            id_part(1) = id_electron
            id_part(2) = id_positron

         ELSEIF (     interaction .EQ. 16
     +           .OR. interaction .EQ. 64) THEN
            npart = 2
            id_part(1) = id_positron
            id_part(2) = id_electron

         ENDIF
      ENDIF

      RETURN
      END
