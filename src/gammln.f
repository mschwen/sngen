      REAL*8 FUNCTION gammln(xx)

      IMPLICIT none

      REAL*8 xx
      INTEGER j
      DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)

      SAVE cof, stp
      DATA cof, stp /76.18009172947146d0, -8650532032941677d0,
     +               24.01409824083091d0, -1.231739572450155d0,
     +               0.1208650973866179d0, -0.5395239384953d-5,
     +               2.5066282746310005d0/

***********************************************************************

       x = xx
       y = x
       tmp = x + 5.5d0
       tmp = (x + 0.5d0) * log(tmp) - tmp
       ser = 1.000000000190015d0

       do j = 1, 6
          y = y + 1.d0
          ser = ser + cof(j) / y
       end do

       gammln = tmp + log(stp * ser / x)

       return
       END
