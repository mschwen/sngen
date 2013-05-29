***********************************************************************
*  Filename: get_int_type.f
***********************************************************************
*
*  This routine returns the interaction class type given a specific
*  interaction:
*  (1:CC, 2:NC, 3:ES)
*  The correspondance between the integer id and character names is done
*  through "id_int.inc".
*   
*  Inputs:
*       o  interaction id number (not array)
*  Outputs:
*       o  neutrino interaction type (integer)
*
***********************************************************************
     
      SUBROUTINE get_int_type(interaction,int_type)

      IMPLICIT none

      INCLUDE 'id_int.inc'
      INCLUDE 'sn_classes.inc'

C  Global Variable Declarations

      INTEGER interaction
      INTEGER int_type

C  Local Variable Declarations

      INTEGER found, i

************************************************************************

      found = 0 

      DO 10 i = 1, cc_size_SNO
         IF (interaction .EQ. list_cc_SNO(i)) THEN
            int_type = id_nu_cc
            found = 1
         ENDIF
10    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 20 i = 1, nc_size_SNO
         IF (interaction .EQ. list_nc_SNO(i)) THEN
            int_type = id_nu_nc
            found = 1
         ENDIF
20    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 30 i = 1, es_size_SNO
         IF (interaction .EQ. list_es_SNO(i)) THEN
            int_type = id_nu_es
            found = 1
         ENDIF
30    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 40 i = 1, cc_size_SNOP
         IF (interaction .EQ. list_cc_SNOP(i)) THEN
            int_type = id_nu_cc
            found = 1
         ENDIF
40    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 50 i = 1, nc_size_SNOP
         IF (interaction .EQ. list_nc_SNOP(i)) THEN
            int_type = id_nu_nc
            found = 1
         ENDIF
50    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 60 i = 1, es_size_SNOP
         IF (interaction .EQ. list_es_SNOP(i)) THEN
            int_type = id_nu_es
            found = 1
         ENDIF
60    CONTINUE

      DO 65 i = 1, es_size_HALO
         IF (interaction .EQ. list_es_HALO(i)) THEN
            int_type = id_nu_es
            found = 1
         ENDIF
65    CONTINUE

      IF (found .EQ. 1) RETURN   

      DO 70 i = 1, cc_size_HALO
         IF (interaction .EQ. list_cc_HALO(i)) THEN
            int_type = id_nu_cc
            found = 1
         ENDIF
70    CONTINUE

      IF (found .EQ. 1) RETURN

      DO 80 i = 1, nc_size_HALO
         IF (interaction .EQ. list_nc_HALO(i)) int_type = id_nu_nc
80    CONTINUE

C      WRITE (*,*) int_type

      RETURN
      END
