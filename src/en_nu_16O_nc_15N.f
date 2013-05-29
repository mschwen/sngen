**************************************************************************
* Filename: en_nu_16O_nc_15N.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_16O_nc_15N(E_part)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION Ran
      DOUBLE PRECISION ran_val
      REAL*8 E_part

*****************************************************************************

      Ran = ran_val(0.d0, 1.d0) 

         IF (Ran .LT. 0.457) THEN
            E_part = 0
         ELSEIF (Ran .LT. 0.523) THEN
            E_part = 5.3
         ELSEIF (Ran .LT. 0.819) THEN
            E_part = 6.33
         ELSEIF (Ran .LT. 0.892) THEN
            E_part = 7.3
         ELSEIF (Ran .LT. 0.919) THEN
            E_part = 9.1
         ELSEIF (Ran .LT. 0.939) THEN
            E_part = 9.2
         ELSEIF (Ran .LT. 0.970) THEN
            E_part = 9.9
         ELSE
            E_part = 10.8
    
         ENDIF

c       WRITE (*,*) E_part

      RETURN
      END
