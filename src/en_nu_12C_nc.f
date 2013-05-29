**************************************************************************
*   Filename: en_nu_12C_nc.f
**************************************************************************
*
*   Outputs:
*       o  sampled final state particle total electron energy (MeV)
*
***************************************************************************

      SUBROUTINE en_nu_12C_nc(E)
   
      IMPLICIT none

C   Local Variable Declarations

      DOUBLE PRECISION Ran
      DOUBLE PRECISION ran_val
      REAL*8 E
      INTEGER n_gamma

      DATA n_gamma /1/
      DATA Ran /0/

      SAVE n_gamma
      SAVE Ran

*****************************************************************************

      IF (Ran .EQ. 0) THEN

            Ran = ran_val(0.d0, 1.d0)

      ENDIF

      IF (n_gamma .EQ. 1) THEN
            IF (Ran .LT. 0.936330) THEN
                E = 15.1
            ELSEIF (Ran .LT. 0.962547) THEN
                E = 7.45
            ELSEIF (Ran .LT. 0.985955) THEN
                E = 10.7 
            ELSEIF (Ran .LT. 0.998270) THEN
                E = 2.40
            ELSE
                E = 2.40
            ENDIF

            n_gamma = 2

      ELSEIF (n_gamma .EQ. 2) THEN
            IF (Ran .LT. 0.936330) THEN
                E = 0 
            ELSEIF (Ran .LT. 0.962547) THEN
                E = 3.21
            ELSEIF (Ran .LT. 0.985955) THEN
                E = 4.44
            ELSEIF (Ran .LT. 0.998270) THEN
                E = 12.7
            ELSE
                E = 8.27
            ENDIF

            n_gamma = 3

      ELSEIF (n_gamma .EQ. 3) THEN
            IF (Ran .LT. 0.936330) THEN
                E = 0 
            ELSEIF (Ran .LT. 0.962547) THEN
                E = 4.44
            ELSEIF (Ran .LT. 0.985955) THEN
                E = 0
            ELSEIF (Ran .LT. 0.998270) THEN
                E = 0
            ELSE
                E = 4.44
            ENDIF

            n_gamma = 1
            Ran = 0

       ENDIF 

       RETURN
       END
