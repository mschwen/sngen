***************************************************************************
*  Filename: poidev.f
***************************************************************************
*
*  Poisson distribution function (taken from Numercial Recipes v2)
*  (replaced "ran1" routine with ral_val routine that calls "ran3";
*   for some reason ran1 was producing errors)
*
***************************************************************************

      FUNCTION poidev(xm,idum)

      IMPLICIT none

      INCLUDE 'trig_constants.inc'

      INTEGER idum
      REAL*8 poidev, xm

C     uses gammln, ran1 (not really)

      REAL*8 alxm, em, g, oldm, sq, t, y, gammln, ran_val
      SAVE alxm, g, oldm, sq
      DATA oldm /-1./

***************************************************************************

      if (xm .lt. 12.) then
         if (xm .ne. oldm) then
            oldm = xm
            g = dexp(-xm)
         endif

         em = -1
         t = 1.d0
2        em = em + 1.
         t = t * ran_val(0.d0,1.d0)
         if (t .gt. g) goto 2
      else
         if (xm .ne. oldm) then
            oldm = xm
            sq = dsqrt(2.d0 * xm)
            alxm = dlog(xm)
            g = xm * alxm - gammln(xm + 1.d0)
         endif

1        y = dtan(pi * ran_val(0.d0,1.d0))
         em = sq * y + xm
         if (em .lt. 0.d0) goto 1
         em = int(em)
         t = 0.9d0 * (1.d0 + y**2) 
     +      * dexp(em * alxm - gammln(em + 1.d0) - g)
         if (ran_val(0.d0,1.d0) .gt. t) goto 1
      endif

      poidev = em

      return
      END
