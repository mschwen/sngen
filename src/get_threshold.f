***********************************************************************
*  Filename: get_threshold.f
***********************************************************************
*
*  This routine specifies the energy threshold for a given reaction
*     
*  Inputs:
*       o  interaction id number (not array)
*  Outputs:
*       o  energy threshold for the specified reaction
*
***********************************************************************
     
      SUBROUTINE get_threshold(interaction,E_thres)

      IMPLICIT none

      INCLUDE 'threshold.inc'

C  Global Variable Declarations

      INTEGER interaction
      DOUBLE PRECISION E_thres

************************************************************************
*  Initialize
************************************************************************

      E_thres = 0.d0

************************************************************************
*  Assign neutrino energy thresholds (MeV)
************************************************************************

      IF (     interaction .EQ. 1
     +    .OR. interaction .EQ. 2
     +    .OR. interaction .EQ. 10
     +    .OR. interaction .EQ. 11
     +    .OR. interaction .EQ. 21
     +    .OR. interaction .EQ. 22
     +    .OR. interaction .EQ. 30
     +    .OR. interaction .EQ. 31
     +    .OR. interaction .EQ. 39
     +    .OR. interaction .EQ. 40
     +    .OR. interaction .EQ. 48
     +    .OR. interaction .EQ. 49
     +    .OR. interaction .EQ. 57
     +    .OR. interaction .EQ. 61
     +    .OR. interaction .EQ. 66
     +    .OR. interaction .EQ. 69
     +    .OR. interaction .EQ. 72
     +    .OR. interaction .EQ. 75
     +    .OR. interaction .EQ. 78
     +    .OR. interaction .EQ. 83
     +    .OR. interaction .EQ. 88
     +    .OR. interaction .EQ. 91
     +    .OR. interaction .EQ. 94
     +    .OR. interaction .EQ. 97) THEN

        E_thres = thres_nu_es_e

      ELSEIF (  interaction .EQ. 58
     +    .OR.  interaction .EQ. 62
     +    .OR.  interaction .EQ. 67
     +    .OR.  interaction .EQ. 70
     +    .OR.  interaction .EQ. 73
     +    .OR.  interaction .EQ. 76) THEN

        E_thres = thres_nu_es_p

      ELSEIF (interaction .EQ. 3) THEN
        E_thres = thres_nu_e_d_cc

      ELSEIF (     interaction .EQ. 4
     +        .OR. interaction .EQ. 59) THEN
        E_thres = thres_nu_e_12C_cc

      ELSEIF (     interaction .EQ. 5
     +        .OR. interaction .EQ. 6
     +        .OR. interaction .EQ. 7) THEN
        E_thres = thres_nu_e_16O_cc

      ELSEIF (     interaction .EQ. 8
     +        .OR. interaction .EQ. 19
     +        .OR. interaction .EQ. 23
     +        .OR. interaction .EQ. 32
     +        .OR. interaction .EQ. 41
     +        .OR. interaction .EQ. 50) THEN
        E_thres = thres_nu_d_nc

      ELSEIF (     interaction .EQ. 9
     +        .OR. interaction .EQ. 20
     +        .OR. interaction .EQ. 24
     +        .OR. interaction .EQ. 33
     +        .OR. interaction .EQ. 42
     +        .OR. interaction .EQ. 51
     +        .OR. interaction .EQ. 60
     +        .OR. interaction .EQ. 65
     +        .OR. interaction .EQ. 68
     +        .OR. interaction .EQ. 71
     +        .OR. interaction .EQ. 74
     +        .OR. interaction .EQ. 77) THEN
        E_thres = thres_nu_12C_nc

      ELSEIF (     interaction .EQ. 12
     +        .OR. interaction .EQ. 13
     +        .OR. interaction .EQ. 14
     +        .OR. interaction .EQ. 63) THEN
        E_thres = thres_nu_e_bar_p_cc

      ELSEIF (interaction .EQ. 15) THEN
        E_thres = thres_nu_e_bar_d_cc

      ELSEIF (     interaction .EQ. 16
     +        .OR. interaction .EQ. 64) THEN
        E_thres = thres_nu_e_bar_12C_cc

      ELSEIF (     interaction .EQ. 17
     +        .OR. interaction .EQ. 18) THEN
        E_thres = thres_nu_e_bar_16O_cc

      ELSEIF (     interaction .EQ. 79) THEN
        E_thres = thres_nu_208Pb_cc_1n

      ELSEIF (     interaction .EQ. 80) THEN
        E_thres = thres_nu_208Pb_cc_2n

      ELSEIF (     interaction .EQ. 81
     +        .OR. interaction .EQ. 89
     +        .OR. interaction .EQ. 95) THEN
       E_thres = thres_nu_208Pb_nc_1n

      ELSEIF (     interaction .EQ. 86
     +        .OR. interaction .EQ. 92
     +        .OR. interaction .EQ. 98) THEN
       E_thres = thres_nub_208Pb_nc_1n

      ELSEIF (     interaction .EQ. 82
     +        .OR. interaction .EQ. 90
     +        .OR. interaction .EQ. 96) THEN
       E_thres = thres_nu_208Pb_nc_2n

      ELSEIF (     interaction .EQ. 87
     +        .OR. interaction .EQ. 93
     +        .OR. interaction .EQ. 99) THEN
       E_thres = thres_nub_208Pb_nc_2n

      ELSE
        E_thres = thres_nu_16O_nc
      ENDIF

c        WRITE(*,*) E_thres

      RETURN
      END
