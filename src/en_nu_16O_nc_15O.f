**************************************************************************
* Filename: en_nu_16O_nc_15O.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_16O_nc_15O(E_part)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION Ran
      DOUBLE PRECISION ran_val
      REAL*8 E_part
      INTEGER n_gamma

      DATA n_gamma /1/

      SAVE n_gamma

*****************************************************************************
      
      IF (n_gamma .EQ. 1) THEN

        E_part = 2.2 ! neutron is captured resulting in 2.2 MeV photon
        n_gamma = 2

      ELSEIF (n_gamma .EQ. 2) THEN

         Ran = ran_val(0.d0, 1.d0) 

         IF (Ran .LT. 0.605) THEN
            E_part = 0
         ELSEIF (Ran .LT. 0.684) THEN
            E_part = 5.2
         ELSEIF (Ran .LT. 0.898) THEN
            E_part = 6.18
         ELSE
            E_part = 6.79
         ENDIF

         n_gamma = 1

      ENDIF

c      WRITE(*,*) E_part

      RETURN
      END
