****************************************************************************
*  Filename: fd_norm.f
****************************************************************************
*
*  This program numerically integrates the Fermi-Dirac function using
*  Simpson's rule.
*
****************************************************************************

	FUNCTION fd_norm(T,eta)

        IMPLICIT none

C     Global Variable Declarations

        REAL*8           T, eta

C     Local Variable Declarations

        REAL*8           fd_norm   ! function
	DOUBLE PRECISION f, x

        DOUBLE PRECISION limit_lo, limit_hi  ! integration limits
        PARAMETER        (limit_lo = 0.d+0)
        PARAMETER        (limit_hi = 5.d+01) ! optimized to agree with
                                             ! Mathematica output
        DOUBLE PRECISION dlimit

        INTEGER          n
        PARAMETER        (n = 5000)           ! must be an even number

        DOUBLE PRECISION dx
        DOUBLE PRECISION xj, odd_sum, even_sum
	INTEGER          j

*****************************************************************************

*****************************************************************************
*  Fermi-Dirac function (needs to be before IF statement)
*****************************************************************************

* changing variables to "x" --> where x = (E_nu / Temp_nu)

        f(x) = (x**2 * T**3) / (1.d0 + dexp(x - eta))

* check to see if integration is necessary

        IF (eta .EQ. 0.d0) THEN   ! the integral is of a standard form
           fd_norm = 1.803085354739584d0 * T**3  ! constant from Mathematica
        ELSE  ! do the integration here

*****************************************************************************
*  Simpson's rule algorithm
*****************************************************************************

           dlimit = (limit_hi - limit_lo)
	   dx = dlimit / dble(n)

	   odd_sum = 0.0
	   DO 10 j = 1, n-1, 2
	      xj = limit_lo + (dx * dble(j))
	      odd_sum = odd_sum + f(xj)
10         CONTINUE

	   even_sum = 0.0
	   DO 20 j = 2, n-2, 2
	      xj = limit_lo + (dx * dble(j))
	      even_sum = even_sum + f(xj)
20         CONTINUE

	   fd_norm = dx / 3.d0 * (f(limit_lo) + f(limit_hi) + 
     +               4.d0 * odd_sum + 2.d0 * even_sum)

        ENDIF

	END
