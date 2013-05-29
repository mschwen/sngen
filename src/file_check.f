**********************************************************************
*  Filename: file_check.f
**********************************************************************
*                                                                     
*  Checks to see if a file exists, returns 0 if not, 1 if yes
*
*  Inputs:
*       o  flag specifying whether or not to open the file if it isn't
*	o  unit_number for OPEN statement
*	o  name of the file to check
*	o  flag indicating whether file exists or not
*
**********************************************************************

      SUBROUTINE file_check(filename,flag)

      IMPLICIT none

C  Global Declarations

      CHARACTER*40 filename
      INTEGER flag

C  Local Declarations

      LOGICAL there

************************************************************************

      flag = 0

      INQUIRE (file = filename, exist = there)

      IF (there) THEN
         flag = 1
      ENDIF

      RETURN
      END
