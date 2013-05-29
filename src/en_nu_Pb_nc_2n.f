* Filename: en_nu_Pb_nc_2n.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_Pb_nc_2n(interaction,E_nu,E_part)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION Ran, ran_val, E_thres
      REAL*8 E_part, E_nu, E_neutron1
      INTEGER n_neutron, interaction

      DATA n_neutron /1/
      DATA E_neutron1 /0/

      SAVE E_neutron1
      SAVE n_neutron

*****************************************************************************

      CALL get_threshold(interaction,E_thres)

      IF (n_neutron .EQ. 1) THEN

          Ran = ran_val(0.d0, 1.d0)

          E_part = (E_nu - E_thres)*Ran

          n_neutron = 2
          E_neutron1 = E_part

      ELSEIF (n_neutron .EQ. 2) THEN

          E_part = E_nu - E_thres - E_neutron1

          n_neutron = 1
          E_neutron1 = 0

      ENDIF

      RETURN
      END
