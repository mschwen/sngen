***********************************************************************
*  Filename: get_xsect.f
***********************************************************************
*
*  Get interaction cross section [m^2]
*
***********************************************************************

      SUBROUTINE get_xsect(interaction,E_nu,xsect_total)
   
      IMPLICIT none

      INCLUDE 'id_int.inc'

C   Global Variable Declarations

      INTEGER           interaction            ! input
      REAL*8            E_nu                   ! input 
      REAL*8            xsect_total            ! output

C   Local Variable Declarations

      INTEGER           int_type

*************************************************************************

      CALL get_int_type(interaction,int_type)

      IF (int_type .EQ. id_nu_es) THEN

         IF (     interaction .EQ. 58 
     +       .OR. interaction .EQ. 62
     +       .OR. interaction .EQ. 67
     +       .OR. interaction .EQ. 70
     +       .OR. interaction .EQ. 73
     +       .OR. interaction .EQ. 76) THEN

              
             CALL xsect_nu_p_es(interaction,E_nu,xsect_total)

         ELSE

             CALL xsect_nu_e_es(interaction,E_nu,xsect_total)

         ENDIF

      ELSEIF (int_type .EQ. id_nu_cc) THEN

********************************************************************
*   3  nu_e + d -> p + p + e-      D2O CC
*  15  nu_e_bar + d -> n + n + e+  D2O CC
********************************************************************

         IF (     interaction .EQ. 3
     +       .OR. interaction .EQ. 15) THEN

            CALL xsect_nu_d_cc(interaction,E_nu,xsect_total)

***********************************************************************
*   4  nu_e + 12C -> 12N + e-      AV CC
*  16  nu_e_bar + 12C -> 12B + e+  AV CC
*  59  nu_e + 12C -> 12N + e-      LAB CC
*  64  nu_e_bar + 12C -> 12B + e+  LAB CC  
***********************************************************************

         ELSEIF (     interaction .EQ. 4
     +           .OR. interaction .EQ. 16
     +           .OR. interaction .EQ. 59
     +           .OR. interaction .EQ. 64) THEN

            CALL xsect_nu_12C_cc(interaction,E_nu,xsect_total)

***********************************************************************
*   5  nu_e + 16O -> 16F + e-      H2O CC
*   6  nu_e + 16O -> 16F + e-      D2O CC 
*   7  nu_e + 16O -> 16F + e-       AV CC 
*  17  nu_e_bar + 16O -> 16N + e+  H2O CC
*  18  nu_e_bar + 16O -> 16N + e+  D2O CC
***********************************************************************

         ELSEIF (     interaction .EQ. 5
     +           .OR. interaction .EQ. 6
     +           .OR. interaction .EQ. 7
     +           .OR. interaction .EQ. 17
     +           .OR. interaction .EQ. 18) THEN

            CALL xsect_nu_16O_cc(interaction,E_nu,xsect_total)

************************************************************************
*  12  nu_e_bar + p -> n + e+     H2O CC
*  13  nu_e_bar + p -> n + e+     D2O* CC  --> H2O in D2O region
*  14  nu_e_bar + p -> n + e+     AV CC
*  63  nu_e_bar + p -> n + e+     LAB CC
************************************************************************

         ELSEIF (     interaction .EQ. 12
     +           .OR. interaction .EQ. 13
     +           .OR. interaction .EQ. 14
     +           .OR. interaction .EQ. 63) THEN

            CALL xsect_nu_p_cc(interaction,E_nu,xsect_total)

************************************************************************
*  79  nu_e + 208Pb -> + 207Bi + n + e-
************************************************************************

         ELSEIF (interaction .EQ. 79) THEN

            CALL xsect_nu_208Pb_cc_1n(interaction,E_nu,xsect_total)

************************************************************************
*  80 nu_e + 208Pb -> 207Bi + 2n + e-
************************************************************************

         ELSEIF (interaction .EQ. 80) THEN

            CALL xsect_nu_208Pb_cc_2n(interaction,E_nu,xsect_total)     

         ELSEIF (     interaction .EQ. 84
     +           .OR. interaction .EQ. 85) THEN

             xsect_total = 0.0

         ENDIF

********
*  Neutral Current
********

      ELSEIF (int_type .EQ. id_nu_nc) THEN

********************************************************************
*   8  nu_e + d --> nu_e + p + n  D2O NC 
*  19  nu_e_bar + d --> nu_e_bar + p + n  D2O NC
*  23  nu_mu + d --> nu_mu + p + n  D2O NC
*  32  nu_mu_bar + d --> nu_mu_bar + p + n  D2O NC
*  41  nu_tau + d --> nu_tau + p + n  D2O NC
*  50  nu_tau_bar + d --> nu_tau_bar + p + n  D2O NC
********************************************************************
      
         IF (     interaction .EQ. 8
     +       .OR. interaction .EQ. 19
     +       .OR. interaction .EQ. 23
     +       .OR. interaction .EQ. 32
     +       .OR. interaction .EQ. 41
     +       .OR. interaction .EQ. 50) THEN

            CALL xsect_nu_d_nc(interaction,E_nu,xsect_total)

********************************************************************
*  These cross section ratios are adapted from the work of 
*   Langanke, Vogel & Kolbe PRL 76 2629 (1996)
*
*  The form for the total cross section is taken from 
*   Beacom & Vogel hep-ph/9802424 v2 10 June 1998
********************************************************************

********************************************************************
*  25  nu_mu + 16O --> nu_mu + 15O + n  D2O NC
*  26  nu_mu + 16O --> nu_mu + 15N + p + gamma  H2O NC
*  27  nu_mu + 16O --> nu_mu + 15N + p + gamma  D2O NC
*  28  nu_mu + 16O --> nu_mu + 15O + n + gamma  H2O NC
*  29  nu_mu + 16O --> nu_mu + 15O + n + gamma  D2O NC
*  34  nu_mu_bar + 16O --> nu_mu_bar + 15O + n  D2O NC
*  35  nu_mu_bar + 16O --> nu_mu_bar + 15N + p + gamma  H2O NC
*  36  nu_mu_bar + 16O --> nu_mu_bar + 15N + p + gamma  D2O NC
*  37  nu_mu_bar + 16O --> nu_mu_bar + 15O + n + gamma  H2O NC
*  38  nu_mu_bar + 16O --> nu_mu_bar + 15O + n + gamma  D2O NC
*  43  nu_tau + 16O --> nu_tau + 15O + n  D2O NC
*  44  nu_tau + 16O --> nu_tau + 15N + p + gamma  H2O NC
*  45  nu_tau + 16O --> nu_tau + 15N + p + gamma  D2O NC
*  46  nu_tau + 16O --> nu_tau + 15O + n + gamma  H2O NC
*  47  nu_tau + 16O --> nu_tau + 15O + n + gamma  D2O NC
*  52  nu_tau_bar + 16O --> nu_tau_bar + 15O + n  D2O NC
*  53  nu_tau_bar + 16O --> nu_tau_bar + 15N + p + gamma  H2O NC
*  54  nu_tau_bar + 16O --> nu_tau_bar + 15N + p + gamma  D2O NC
*  55  nu_tau_bar + 16O --> nu_tau_bar + 15O + n + gamma  H2O NC
*  56  nu_tau_bar + 16O --> nu_tau_bar + 15O + n + gamma  D2O NC
********************************************************************

         ELSEIF (     interaction .EQ. 25
     +           .OR. interaction .EQ. 26
     +           .OR. interaction .EQ. 27
     +           .OR. interaction .EQ. 28
     +           .OR. interaction .EQ. 29
     +           .OR. interaction .EQ. 34
     +           .OR. interaction .EQ. 35
     +           .OR. interaction .EQ. 36
     +           .OR. interaction .EQ. 37
     +           .OR. interaction .EQ. 38
     +           .OR. interaction .EQ. 43
     +           .OR. interaction .EQ. 44
     +           .OR. interaction .EQ. 45
     +           .OR. interaction .EQ. 46
     +           .OR. interaction .EQ. 47
     +           .OR. interaction .EQ. 52
     +           .OR. interaction .EQ. 53
     +           .OR. interaction .EQ. 54
     +           .OR. interaction .EQ. 55
     +           .OR. interaction .EQ. 56) THEN

            CALL xsect_nu_16O_nc(interaction,E_nu,xsect_total)

********************************************************************
*   9  nu_e + 12C -> nu_e + 12C* (-> gamma)  AV NC
*  20  nu_e_bar + 12C -> nu_e_bar + 12C* (-> gamma)  AV NC
*  24  nu_mu + 12C -> nu_mu + 12C* (-> gamma)  AV NC
*  33  nu_mu_bar + 12C -> nu_mu_bar + 12C* (-> gamma)  AV NC
*  42  nu_tau + 12C -> nu_tau + 12C* (-> gamma)  AV NC
*  51  nu_tau_bar + 12C -> nu_tau_bar + 12C* (-> gamma)  AV NC
*  60  nu_e + 12C -> nu_e + 12C* (-> gamma)  LAB NC
*  65  nu_e_bar + 12C -> nu_e_bar + 12C* (-> gamma) LAB NC
*  68  nu_mu + 12C -> nu_mu + 12C* (-> gamma)  LAB NC
*  71  nu_mu_bar + 12C -> nu_mu_bar + 12C* (-> gamma)  LAB NC
*  74  nu_tau + 12C -> nu_tau + 12C* (-> gamma)  LAB NC
*  77  nu_tau_bar + 12C -> nu_tau_bar + 12C* (-> gamma)  LAB NC
********************************************************************
      
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
         
            CALL xsect_nu_12C_nc(interaction,E_nu,xsect_total)

*********************************************************************
*  81  nu_e + 208Pb -> nu_e + 207Pb +  n  Pb NC
*  86  nu_e_bar + 208Pb -> nu_ebar + 207Pb +  n  Pb NC
*  89  nu_mu + 208Pb -> nu_mu + 207Pb + n  Pb NC
*  92  nu_mu_bar + 208Pb -> nu_mu_bar + 207Pb + n  Pb NC
*  95  nu_tau + 208Pb -> nu_tau + 207Pb + n  Pb NC
*  98  nu_tau_bar + 208Pb -> nu_tau_bar + 207Pb + n  Pb NC
*********************************************************************

        ELSEIF (     interaction .EQ. 81
     +          .OR. interaction .EQ. 86
     +          .OR. interaction .EQ. 89
     +          .OR. interaction .EQ. 92
     +          .OR. interaction .EQ. 95
     +          .OR. interaction .EQ. 98) THEN
   
            CALL xsect_nu_208Pb_nc_1n(interaction,E_nu,xsect_total)

*********************************************************************
*  82  nu_e + 208Pb -> nu_e + 207Pb + 2n  Pb NC
*  87  nu_e_bar + 208Pb -> nu_ebar + 207Pb +  2n  Pb NC
*  90  nu_mu + 208Pb -> nu_mu + 207Pb + 2n  Pb NC
*  93  nu_mu_bar + 208Pb -> nu_mu_bar + 207Pb + 2n  Pb NC
*  96  nu_tau + 208Pb -> nu_tau + 207Pb + 2n  Pb NC
*  99  nu_tau_bar + 208Pb -> nu_tau_bar + 207Pb + 2n  Pb NC
*********************************************************************

        ELSEIF (     interaction .EQ. 82
     +          .OR. interaction .EQ. 87
     +          .OR. interaction .EQ. 90
     +          .OR. interaction .EQ. 93
     +          .OR. interaction .EQ. 96
     +          .OR. interaction .EQ. 99) THEN

            CALL xsect_nu_208Pb_nc_2n(interaction,E_nu,xsect_total)


         ENDIF  ! end of neutral current
      ENDIF     ! end of interaction classes

C      WRITE (*,*) E_nu, xsect_total

      RETURN
      END

