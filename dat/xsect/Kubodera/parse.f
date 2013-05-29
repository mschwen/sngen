      PROGRAM parse_Kubodera

      INTEGER nentries, i, j
      PARAMETER (nentries = 136)

      REAL*8  Enu, xs_num(4)

      CHARACTER*5 xs_char(4)

***************************************************************************

      OPEN (unit = 10, file = 'nsapmgk.dat', status = 'old')
      OPEN (unit = 20, file = 'xs12.dat', status = 'unknown')
      OPEN (unit = 30, file = 'xs34.dat', status = 'unknown')

      DO i = 1, nentries
         READ (10,*) Enu

         DO j = 1, 4
            READ (10,*) xs_num(j), xs_char(j)
         END DO

         IF (Enu .EQ. 0.) THEN
            WRITE (20,10) Enu, xs_num(1), xs_num(2)
            WRITE (30,10) Enu, xs_num(3), xs_num(4)
         ELSE
            WRITE (20,20) Enu, xs_num(1), xs_char(1), 
     +                    xs_num(2), xs_char(2)
            WRITE (30,20) Enu, xs_num(3), xs_char(3), 
     +                    xs_num(4), xs_char(4)
         ENDIF 
      END DO

10    FORMAT (f6.2, 3x, f5.3, 3x, f5.3)
20    FORMAT (f6.2, 3x, f5.3,'e', a, 3x, f5.3,'e', a)

      CLOSE (unit = 20)
      CLOSE (unit = 30)

      STOP
      END
