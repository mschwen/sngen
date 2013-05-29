************************************************************************
*  Filename: aspline_eval.f
************************************************************************
*
*  This routine is used to evaluate the coefficients calculated in
*  the 'aspline' routine.  This is a modified version of Paul Turner's 
*  'seval' routine in "fit.c" that is part of the plotting package "xmgr". 
*
*   Inputs:
*	o  ndata = initial number of data points to spline fit
*	o  given 'x' variable data array
*	o  given 'y' variable data array
*	o  interpolation range in x (xstart,xstop)
*	o  number of spline interpolation points
*
************************************************************************

      SUBROUTINE aspline_eval(n_data,x,y,b,c,d,xfit,yfit)

      IMPLICIT none

C   Global Variable Declarations

      INTEGER n_data                                   ! input
      DOUBLE PRECISION x(*), y(*)                      ! input
      DOUBLE PRECISION b(*), c(*), d(*)                ! input
      DOUBLE PRECISION xfit, yfit                      ! output

C   Local Variable Declarations

      DOUBLE PRECISION dx, dy, m                       ! for n = 1, 2 case
      INTEGER khi, klo, k

**************************************************************************
*   Initializations
**************************************************************************

      khi = 0
      klo = 0

**************************************************************************
*  Treat the n_data = 1, n_data = 2 as special cases
**************************************************************************

      IF (n_data .GT. 0) THEN

         IF (n_data .EQ. 1) THEN
            xfit = x(1)
            yfit = y(1)
          
         ELSEIF (n_data .EQ. 2) THEN
            dx = x(2) - x(1)
            dy = y(2) - y(1)

            m = dy / dx

            yfit = m * xfit

***************************************************************************
*  We have more than two data points with which to evaluate the spline fit
***************************************************************************

         ELSE

* binary search with specified 'xfit' value

            IF (xfit .LT. x(1)) THEN
               klo = 1
            ELSEIF (xfit .GE. x(n_data)) THEN
               klo = n_data - 1
            ELSE
               klo = 1
               khi = n_data

               DOWHILE ((khi - klo) .GT. 1)
                  k = (khi + klo) / 2
                  IF (xfit .LT. x(k)) THEN
                     khi = k
                  ELSE
                     klo = k
                  ENDIF
               END DO
            ENDIF

            dx = xfit - x(klo)

*  calculate Akima spline fitted value at the specified 'xfit' value

            yfit = y(klo) + dx * (b(klo) + dx * (c(klo) + dx * d(klo)))

         ENDIF
      ENDIF

      RETURN
      END
