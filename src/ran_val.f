**********************************************************************
*  Filename: ran_val.f
**********************************************************************
*
*   Generates random value between "max" and "min"
*
*   Inputs:
*	o  minimum variable value
*	o  maximum variable value
*
**********************************************************************
        
      DOUBLE PRECISION FUNCTION ran_val(min,max)

      IMPLICIT none

      INCLUDE 'random_seed.inc'

C  Global Variable Declarations

      DOUBLE PRECISION min, max

C  Local Variable Declarations

      DOUBLE PRECISION ran, ran3
      DOUBLE PRECISION diff

      INTEGER count

      SAVE count
      DATA count /0/

***********************************************************************

      ran = 0.

***********************************************************************
*  Choose random value between max and min.
*  The first call to 'ran3' initializes an array of random numbers
*  generated with 'iseed' as the starting value.
*  Subsequent calls require '1' as the argument to access successive
*  random numbers in the array.
***********************************************************************

      IF (count .EQ. 0) THEN
         ran = ran3(iseed)
         count = 1
      ELSE
         ran = ran3(1)
      ENDIF

      diff = max - min
      ran_val = min + (diff * ran)

      IF (ran_val .LT. min .OR. ran_val .GT. max) THEN
         WRITE (*,*) '*** ran_val.f: Something bad happened!'
         WRITE (*,*) '*** ran_val.f: ran_val =', ran_val
      ENDIF

      RETURN
      END
