*FILE MEMBER=MZCA LIBRARY=SNOMAN LANGUAGE=INC DATE=26:May:1997
 
*     Standard ZEBRA Common:  DO NOT MODIFY!
 
*     Contact: R. Komar, UBC.
 
*     Revision History:-
*     ================
 
*     3_00  Rob Komar   First version.
c
c  The following is a Zebra system common block, which we use to check
c  to see if Zebra has already been initialized.
c
      INTEGER NQSTOR,NQOFFT,NQOFFS,NQALLO, NQIAM
     +,              LQATAB,LQASTO,LQBTIS, LQWKTB,NQWKTB,LQWKFZ
     +,              MQKEYS,NQINIT,NQTSYS,NQM99,NQPERM,NQFATA,NQCASE
     +,              NQTRAC,MQTRAC
      COMMON /MZCA/  NQSTOR,NQOFFT(16),NQOFFS(16),NQALLO(16), NQIAM
     +,              LQATAB,LQASTO,LQBTIS, LQWKTB,NQWKTB,LQWKFZ
     +,              MQKEYS(3),NQINIT,NQTSYS,NQM99,NQPERM,NQFATA,NQCASE
     +,              NQTRAC,MQTRAC(48)
      INTEGER KQSP
                                       EQUIVALENCE (KQSP,NQOFFS(1))
      INTEGER INKEYS(3)
      DATA  INKEYS / 4HEBRA, 4HINIT, 4HDONE /
 
*ENDFILE MEMBER=MZCA
