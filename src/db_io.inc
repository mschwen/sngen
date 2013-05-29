*FILE MEMBER=DB_IO LIBRARY=SNODB LANGUAGE=INC DATE=23:Apr:1997
 
*     SNODB:  logical unit numbers for SNODB internal use.
 
*     Contact:  R. Komar (UBC)
*               D. Cowen (UPenn)
 
*     Revision History:-
*     ================
 
*     2.09  R. Komar     First version
*     3.00  R. Komar	 Replaced all luncdxx with single luncd
 
 
*     Variables:-
*     =========
 
	integer luncd	!logical unit number for database file
	integer lunfz	!logical unit number for server updates
	integer	luntx	!logical unit number for text files
	integer luntz	!logical unit number for titles files
	integer luntzlist	!logical unit number for temporary files
c
	integer lrecl	!record length for database file
	integer luntmp
 
	parameter (luntz=40)
	parameter (luntx=41)
	parameter (luntzlist=42)
 
*ENDFILE MEMBER=DB_IO
