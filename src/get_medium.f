***********************************************************************
*  Filename: get_medium.f
***********************************************************************
*
*  Given the interaction id number, this routine returns the following
*  information:
*
*  Inputs:
*	o interaction id number (not array)
*
*  Outputs:
*	o medium code needed by SNOMAN to specify event position
*
*********************************************************************** 

      SUBROUTINE get_medium(interaction,medium)

      IMPLICIT none

C  Global Variable Declarations

      INTEGER interaction
      CHARACTER*40 medium

**********************************************************************

      IF (     interaction .EQ. 1
     +    .OR. interaction .EQ. 5
     +    .OR. interaction .EQ. 10 
     +    .OR. interaction .EQ. 12
     +    .OR. interaction .EQ. 17
     +    .OR. interaction .EQ. 21
     +    .OR. interaction .EQ. 26
     +    .OR. interaction .EQ. 28
     +    .OR. interaction .EQ. 30
     +    .OR. interaction .EQ. 35
     +    .OR. interaction .EQ. 37
     +    .OR. interaction .EQ. 39
     +    .OR. interaction .EQ. 44
     +    .OR. interaction .EQ. 46
     +    .OR. interaction .EQ. 48
     +    .OR. interaction .EQ. 53
     +    .OR. interaction .EQ. 55) THEN

         medium = '1. $H2O$0000. 1.'

      ELSEIF (     interaction .EQ. 4
     +        .OR. interaction .EQ. 7
     +        .OR. interaction .EQ. 9
     +        .OR. interaction .EQ. 14
     +        .OR. interaction .EQ. 16
     +        .OR. interaction .EQ. 20
     +        .OR. interaction .EQ. 24
     +        .OR. interaction .EQ. 33
     +        .OR. interaction .EQ. 42
     +        .OR. interaction .EQ. 51) THEN

*** acrylic volume is defined as the acrylic sphere + chimney 
         medium = '2. $ACRV$0000. 1. $ACRC$0000. 1.'

      ELSE

*** D2O volume is define by the acrylic sphere inner volume 
*** + chimney inner volume
         medium = '2. $D2O$0000. 1. $ACRC_IVL$0000. 1.'
      ENDIF

      RETURN
      END
