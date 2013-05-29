***********************************************************************
*  Filename: file_cmd_type.f
***********************************************************************
*
*  Appends a modifier on certain cmd files identifying whether or not
*  the cmd file is for a "salt" or "nosalt" run.
*
*  Note: The default filename is assumed to be of the form: 
*        	<file>##.out
*        Basically, the routine looks for the '.' in the default_name
*        character string and adds the type labels at the position of
*        the '.'.
*
*  Inputs:
*	o  default filename to be modified
*	o  the increment by which the ## will be modified
*
***********************************************************************

      SUBROUTINE file_cmd_type(default_name,type,output_name)

      IMPLICIT none

      INCLUDE 'id_cmd_type.inc'

C  Global Variable Declarations

      INTEGER type
      CHARACTER*(*) default_name
      CHARACTER*40 output_name

C  Local Variable Declarations

      CHARACTER*40 dummy_name
      INTEGER length_part1, length_total

************************************************************************

************************************************************************
*  Initializations
************************************************************************

      length_total = 0
      length_part1 = 0
c      dummy_name = default_name

*************************************************************************

      length_part1 = index(default_name, '.')
      length_total = index(default_name, ' ')

      IF (     length_total .EQ. 0
     +    .OR. length_part1 .EQ. 0) THEN
         WRITE (*,*) 'The following file is not of the correct format:',
     +   default_name
      ENDIF 

c      WRITE (*,*) 'filename, length =', default_name, length_total

      IF (type .EQ. id_nosalt) THEN
         dummy_name = default_name(1:length_part1-1) // "_nosalt" //
     + default_name(length_part1:length_total)

      ELSEIF (type .EQ. id_salt) THEN
         dummy_name = default_name(1:length_part1-1) // "_salt" //
     + default_name(length_part1:length_total)

      ELSEIF (type .EQ. id_ncd) THEN
         dummy_name = default_name(1:length_part1-1) // "_ncd" //
     + default_name(length_part1:length_total)

      ELSE
         WRITE (*,*) '*** file_cmd_type.f: Something went wrong!'
      ENDIF

      output_name = dummy_name

      RETURN
      END
