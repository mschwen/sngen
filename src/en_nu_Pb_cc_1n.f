* Filename: en_nu_Pb_cc_1n.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_Pb_cc_1n(interaction,E_nu,E_part)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION E_thres
      REAL*8 E_part, E_nu
      INTEGER interaction

*****************************************************************************

      CALL get_threshold(interaction,E_thres)

      E_part = (E_nu - E_thres)

      RETURN
      END
