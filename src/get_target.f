***********************************************************************
*  Filename: get_target.f
***********************************************************************
*
*  This routine specifies the energy threshold for a given reaction
*     
*  Inputs:
*       o  interaction id number (not array)
*  Outputs:
*       o  number of target atoms for the specified reaction
*
***********************************************************************
     
      SUBROUTINE get_target(interaction,N_target)

      IMPLICIT none

      INCLUDE 'targets.inc'

C  Common Block Declarations

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

      COMMON           /sn_model_params/
     +                 use_author_params
      INTEGER          use_author_params

C  Global Variable Declarations

      INTEGER          interaction  ! input
      REAL*8           N_target     ! output

C  Local Variable Declarations

************************************************************************

      IF (     interaction .EQ. 1
     +    .OR. interaction .EQ. 10
     +    .OR. interaction .EQ. 21
     +    .OR. interaction .EQ. 30
     +    .OR. interaction .EQ. 39
     +    .OR. interaction .EQ. 48) THEN

        N_target = N_e_h2o

        IF (use_author_params .EQ. 1) THEN
           IF (sn_model_switch .EQ. 1) THEN   ! Burrows
              N_target = (1.6d0 / 1.7d0) * N_target
           ELSEIF (sn_model_switch .EQ. 2) THEN  ! Beacom & Vogel
              N_target = (1.4d0 / 1.7d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (     interaction .EQ. 2
     +        .OR. interaction .EQ. 11
     +        .OR. interaction .EQ. 22
     +        .OR. interaction .EQ. 31
     +        .OR. interaction .EQ. 40
     +        .OR. interaction .EQ. 49) THEN

        N_target = N_e_d2o

        IF (use_author_params .EQ. 1) THEN
           IF (     sn_model_switch .EQ. 1   ! Burrows
     +         .OR. sn_model_switch .EQ. 2)  ! Beacom & Vogel
     +     THEN
              N_target = (30.08d31 / N_e_d2o) * N_target  ! from Burrows
c              N_target = (1.d0 / 0.9992d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (     interaction .EQ. 3
     +        .OR. interaction .EQ. 8
     +        .OR. interaction .EQ. 15
     +        .OR. interaction .EQ. 19
     +        .OR. interaction .EQ. 23
     +        .OR. interaction .EQ. 32
     +        .OR. interaction .EQ. 41
     +        .OR. interaction .EQ. 50) THEN

        N_target = N_deuteron

        IF (use_author_params .EQ. 1) THEN
           IF (     sn_model_switch .EQ. 1   ! Burrows
     +         .OR. sn_model_switch .EQ. 2)  ! Beacom & Vogel
     +     THEN   
              N_target = (6.016d31 / N_deuteron) * N_target ! from Burrows
c              N_target = (1.d0 / 0.9992d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (     interaction .EQ. 4
     +        .OR. interaction .EQ. 9
     +        .OR. interaction .EQ. 16
     +        .OR. interaction .EQ. 20
     +        .OR. interaction .EQ. 24
     +        .OR. interaction .EQ. 33
     +        .OR. interaction .EQ. 42
     +        .OR. interaction .EQ. 51) THEN

        N_target = N_12C_AV

      ELSEIF (     interaction .EQ. 5
     +        .OR. interaction .EQ. 17
     +        .OR. interaction .EQ. 26
     +        .OR. interaction .EQ. 28
     +        .OR. interaction .EQ. 35
     +        .OR. interaction .EQ. 37
     +        .OR. interaction .EQ. 44
     +        .OR. interaction .EQ. 46
     +        .OR. interaction .EQ. 53
     +        .OR. interaction .EQ. 55) THEN

        N_target = N_16O_h2o

        IF (use_author_params .EQ. 1) THEN
           IF (sn_model_switch .EQ. 1) THEN   ! Burrows
              N_target = (1.6d0 / 1.7d0) * N_target
           ELSEIF (sn_model_switch .EQ. 2) THEN  ! Beacom & Vogel
              N_target = (1.4d0 / 1.7d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (     interaction .EQ. 6
     +        .OR. interaction .EQ. 18
     +        .OR. interaction .EQ. 25
     +        .OR. interaction .EQ. 27
     +        .OR. interaction .EQ. 29
     +        .OR. interaction .EQ. 34
     +        .OR. interaction .EQ. 36
     +        .OR. interaction .EQ. 38
     +        .OR. interaction .EQ. 43
     +        .OR. interaction .EQ. 45
     +        .OR. interaction .EQ. 47
     +        .OR. interaction .EQ. 52
     +        .OR. interaction .EQ. 54
     +        .OR. interaction .EQ. 56) THEN

        N_target = N_16O_d2o

        IF (use_author_params .EQ. 1) THEN
           IF (     sn_model_switch .EQ. 1   ! Burrows
     +         .OR. sn_model_switch .EQ. 2)  ! Beacom & Vogel
     +     THEN   
              N_target = (3.008d31 / N_16O_d2o) * N_target  ! from Burrows
c              N_target = (1.d0 / 0.9992d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (interaction .EQ. 7) THEN

        N_target = N_16O_av

      ELSEIF (interaction .EQ. 12) THEN

        N_target = N_proton_h2o

        IF (use_author_params .EQ. 1) THEN
           IF (sn_model_switch .EQ. 1) THEN   ! Burrows
              N_target = (1.6d0 / 1.7d0) * N_target
           ELSEIF (sn_model_switch .EQ. 2) THEN  ! Beacom & Vogel
              N_target = (1.4d0 / 1.7d0) * N_target
           ENDIF
        ENDIF

      ELSEIF (interaction .EQ. 13) THEN

        N_target = N_proton_d2o

        IF (use_author_params .EQ. 1) THEN
           IF (     sn_model_switch .EQ. 1   ! Burrows
     +         .OR. sn_model_switch .EQ. 2)  ! Beacom & Vogel
     +     THEN   
              N_target = (6.016d31 / N_deuteron) * N_target  ! from Burrows
           ENDIF
        ENDIF

      ELSEIF (interaction .EQ. 14) THEN

        N_target = N_proton_av

      ELSEIF (    interaction .EQ. 57
     +       .OR. interaction .EQ. 61
     +       .OR. interaction .EQ. 66
     +       .OR. interaction .EQ. 69
     +       .OR. interaction .EQ. 72
     +       .OR. interaction .EQ. 75) THEN

        N_target = N_e_LAB
   
      ELSEIF (interaction .EQ. 58
     +       .OR. interaction .EQ. 62
     +       .OR. interaction .EQ. 63
     +       .OR. interaction .EQ. 67
     +       .OR. interaction .EQ. 70
     +       .OR. interaction .EQ. 73
     +       .OR. interaction .EQ. 76) THEN

        N_target = N_proton_LAB

      ELSEIF (    interaction .EQ. 59
     +       .OR. interaction .EQ. 60
     +       .OR. interaction .EQ. 64
     +       .OR. interaction .EQ. 65
     +       .OR. interaction .EQ. 68
     +       .OR. interaction .EQ. 71
     +       .OR. interaction .EQ. 74
     +       .OR. interaction .EQ. 77) THEN

       N_target = N_12C_LAB

      ELSEIF (    interaction .EQ. 78
     +       .OR. interaction .EQ. 79
     +       .OR. interaction .EQ. 80
     +       .OR. interaction .EQ. 81
     +       .OR. interaction .EQ. 82
     +       .OR. interaction .EQ. 83
     +       .OR. interaction .EQ. 84
     +       .OR. interaction .EQ. 85
     +       .OR. interaction .EQ. 86
     +       .OR. interaction .EQ. 87
     +       .OR. interaction .EQ. 88
     +       .OR. interaction .EQ. 89
     +       .OR. interaction .EQ. 90
     +       .OR. interaction .EQ. 91
     +       .OR. interaction .EQ. 92
     +       .OR. interaction .EQ. 93
     +       .OR. interaction .EQ. 94
     +       .OR. interaction .EQ. 95
     +       .OR. interaction .EQ. 96
     +       .OR. interaction .EQ. 97
     +       .OR. interaction .EQ. 98
     +       .OR. interaction .EQ. 99) THEN

      N_target = N_Pb

      ENDIF

      RETURN
      END
