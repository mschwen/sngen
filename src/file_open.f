***********************************************************************
*  Filename: file_open.f
***********************************************************************
*                                                                     
*  Increments the file according to the interaction number 
*
*  Note: the default filename is assumed to be of the form: 
*        <file>##.out
*        (basically, the routine looks for the '.' in the default_name
*         character string and replaces the appropriate # with the
*         increment id number)
*
*  Inputs:
*	o  unit number for the OPEN statement
*	o  default filename to be modified
*	o  the increment by which the ## will be modified
*
***********************************************************************

      SUBROUTINE file_open(unit_number,default_name,increment)

      IMPLICIT none

C  Global Variable Declarations

      INTEGER unit_number
      INTEGER increment
      CHARACTER*(*) default_name

C  Local Variable Declarations

      CHARACTER*40 dummy_name
      CHARACTER*40 output_name
      INTEGER length

************************************************************************

************************************************************************
*  Initializations
************************************************************************

      length = 0
      dummy_name = default_name

*************************************************************************

      length = index(dummy_name, '.')

      IF (length .EQ. 0) THEN
         WRITE (*,*) 'The following file is not of the correct format:',
     +   default_name
      ENDIF 

c      WRITE (*,*) 'filename, length =', default_name, length

      IF (increment .LT. 10) THEN
         WRITE (dummy_name(length-1:length-1), '(i1)') increment
      ELSEIF (increment .GE. 10) THEN
         WRITE (dummy_name(length-2:length-1), '(i2)') increment
      ELSE
         WRITE (*,*) '*** file_open.f: Something went wrong!'
      ENDIF

      output_name = dummy_name

      OPEN (unit=unit_number, file=output_name, status='unknown')

      RETURN
      END
