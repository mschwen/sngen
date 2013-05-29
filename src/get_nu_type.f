***********************************************************************
*  Filename: get_nu_type.f
***********************************************************************
*
*  This routine returns the neutrino type given a specific interaction
*  (the actual integer code is the same one that SNOMAN uses).
*  The correspondance between the integer id and character names is done
*  through "id_nu.inc".
*   
*  Inputs:
*       o  interaction id number (not array)
*  Outputs:
*       o  neutrino type (integer)
*
***********************************************************************
     
      SUBROUTINE get_nu_type(interaction,nu_type)

      IMPLICIT none

      INCLUDE 'sn_classes.inc'
      INCLUDE 'id_nu.inc'

C  Common Block Declarations

C      COMMON /nu_class/ list_nu_e_SNO, list_nu_eb_SNO, 
C     +                  list_nu_mu_SNO, list_nu_mub_SNO,
C     +                  list_nu_tau_SNO, list_nu_taub_SNO,
C     +                  list_nu_e_SNOP, list_nu_eb_SNOP,
C     +                  list_nu_mu_SNOP, list_nu_mub_SNOP,
C     +                  list_nu_tau_SNOP, list_nu_taub_SNOP,
C     +                  list_nu_e_HALO, list_nu_eb_HALO,
C     +                  list_nu_mu_HALO, list_nu_mub_HALO,
C     +                  list_nu_tau_HALO, list_nu_taub_HALO

C      INTEGER list_nu_e_SNO(nu_e_size_SNO), 
C     +        list_nu_eb_SNO(nu_eb_size_SNO),
C     +        list_nu_mu_SNO(nu_mu_size_SNO), 
C     +        list_nu_mub_SNO(nu_mub_size_SNO),
C     +        list_nu_tau_SNO(nu_tau_size_SNO),
C     +        list_nu_taub_SNO(nu_taub_size_SNO),
C     +        list_nu_e_SNOP(nu_e_size_SNOP), 
C     +        list_nu_eb_SNOP(nu_eb_size_SNOP),
C     +        list_nu_mu_SNOP(nu_mu_size_SNOP),
C     +        list_nu_mub_SNOP(nu_mub_size_SNOP),
C     +        list_nu_tau_SNOP(nu_tau_size_SNOP), 
C     +        list_nu_taub_SNOP(nu_taub_size_SNOP),
C     +        list_nu_e_HALO(nu_e_size_HALO), 
C     +        list_nu_eb_HALO(nu_eb_size_HALO),
C     +        list_nu_mu_HALO(nu_mu_size_HALO),
C     +        list_nu_mub_HALO(nu_mub_size_HALO),
C     +        list_nu_tau_HALO(nu_tau_size_HALO), 
C     +        list_nu_taub_HALO(nu_taub_size_HALO)


C  Global Variable Declarations

      INTEGER interaction
      INTEGER nu_type

C  Local Variable Declarations

      INTEGER i

************************************************************************

      DO 10 i = 1, nu_e_size_SNO
         IF (interaction .EQ. list_nu_e_SNO(i)) THEN
            nu_type = id_nu_e
            RETURN
         ENDIF
10    CONTINUE


      DO 20 i = 1, nu_eb_size_SNO
         IF (interaction .EQ. list_nu_eb_SNO(i)) THEN
            nu_type = id_nu_ebar
            RETURN
         ENDIF
20    CONTINUE

      DO 30 i = 1, nu_mu_size_SNO
         IF (interaction .EQ. list_nu_mu_SNO(i)) THEN
            nu_type = id_nu_mu
            RETURN
         ENDIF
30    CONTINUE

      DO 40 i = 1, nu_mub_size_SNO
         IF (interaction .EQ. list_nu_mub_SNO(i)) THEN
            nu_type = id_nu_mubar
            RETURN
         ENDIF
40    CONTINUE

      DO 50 i = 1, nu_tau_size_SNO
         IF (interaction .EQ. list_nu_tau_SNO(i)) THEN
            nu_type = id_nu_tau
            RETURN
         ENDIF
50    CONTINUE

      DO 60 i = 1, nu_taub_size_SNO
         IF (interaction .EQ. list_nu_taub_SNO(i)) THEN
            nu_type = id_nu_taubar
            RETURN
         ENDIF
60    CONTINUE

      DO 70 i = 1, nu_e_size_SNOP
         IF (interaction .EQ. list_nu_e_SNOP(i)) THEN
            nu_type = id_nu_e
            RETURN
         ENDIF
70    CONTINUE

      DO 80 i = 1, nu_eb_size_SNOP
         IF (interaction .EQ. list_nu_eb_SNOP(i)) THEN
            nu_type = id_nu_ebar
            RETURN
         ENDIF
80    CONTINUE

      DO 90 i = 1, nu_mu_size_SNOP
         IF (interaction .EQ. list_nu_mu_SNOP(i)) THEN
            nu_type = id_nu_mu
            RETURN
         ENDIF
90    CONTINUE

      DO 100 i = 1, nu_mub_size_SNOP
         IF (interaction .EQ. list_nu_mub_SNOP(i)) THEN
            nu_type = id_nu_mubar
            RETURN
         ENDIF
100    CONTINUE

      DO 110 i = 1, nu_tau_size_SNOP
         IF (interaction .EQ. list_nu_tau_SNOP(i)) THEN
            nu_type = id_nu_tau
            RETURN
         ENDIF
110    CONTINUE

      DO 120 i = 1, nu_taub_size_SNOP
         IF (interaction .EQ. list_nu_taub_SNOP(i)) THEN
            nu_type = id_nu_taubar
            RETURN
         ENDIF
120    CONTINUE

      DO 130 i = 1, nu_e_size_HALO
         IF (interaction .EQ. list_nu_e_HALO(i)) THEN
            nu_type = id_nu_e
            RETURN
         ENDIF
130    CONTINUE

      DO 140 i = 1, nu_eb_size_HALO
         IF (interaction .EQ. list_nu_eb_HALO(i)) THEN
            nu_type = id_nu_ebar
            RETURN
         ENDIF
140    CONTINUE

      DO 150 i = 1, nu_mu_size_HALO
         IF (interaction .EQ. list_nu_mu_HALO(i)) THEN
            nu_type = id_nu_mu
            RETURN
         ENDIF
150    CONTINUE

      DO 160 i = 1, nu_mub_size_HALO
         IF (interaction .EQ. list_nu_mub_HALO(i)) THEN
            nu_type = id_nu_mubar
            RETURN
         ENDIF
160    CONTINUE

      DO 170 i = 1, nu_tau_size_HALO
         IF (interaction .EQ. list_nu_tau_HALO(i)) THEN
            nu_type = id_nu_tau
            RETURN
         ENDIF
170    CONTINUE

      DO 180 i = 1, nu_taub_size_HALO
         IF (interaction .EQ. list_nu_taub_HALO(i)) THEN
            nu_type = id_nu_taubar
            RETURN
         ENDIF
180    CONTINUE


      RETURN
      END
