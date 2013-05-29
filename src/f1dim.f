****************************************************************************
*  Filename: f1dim.f
****************************************************************************
*
*     The following function for minimization is used by the
*     spectrum generator.  It is taken straight from Numerical
*     Recipes.  It is called from both MNBRAK and BRENT2, and uses
*     an EN_NC/EN_CC specific function called F_DIFF.
* 
****************************************************************************

      FUNCTION f1dim(x)

      IMPLICIT none

      INTEGER NMAX
      REAL f1dim,xsect_nu_d_diff,x
      PARAMETER (NMAX = 50)
CU    USES xsect_nu_d_diff
      INTEGER j, ncom
      REAL pcom(NMAX), xicom(NMAX), xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom

*****************************************************************************

      DO 11 j = 1, ncom
        xt(j) = pcom(j) + x * xicom(j)
11    CONTINUE

      f1dim = xsect_nu_d_diff(xt)

      RETURN
      END
