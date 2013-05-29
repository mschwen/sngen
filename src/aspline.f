************************************************************************
*  Filename: aspline.f
************************************************************************
*
*   Akima spline fitter subroutine 
*   (modified version of David Frey's "spline.c" /
*    Paul Turner's "fit.c" aspline routine used with "xmgr") 
*
*   Inputs:
*	o  n_data = initial number of data points to spline fit
*	o  given 'x' variable data array
*	o  given 'y' variable data array
*
*   Outputs:
*	o  b,c,d are Akima spline coefficients
*
************************************************************************

      SUBROUTINE aspline(n_data,x,y,b,c,d)

      IMPLICIT none

      INCLUDE 'aspline_size.inc'

C   Global Variable Declarations

      INTEGER n_data                                          ! input
      DOUBLE PRECISION x(*), y(*)                             ! input
      DOUBLE PRECISION b(aspline_size), c(aspline_size),
     +                 d(aspline_size)                        !output

C   Local Variable Declarations

      DOUBLE PRECISION dx(aspline_size), dy(aspline_size), 
     +                 m(aspline_size)
      DOUBLE PRECISION num, den
      DOUBLE PRECISION m_m1, m_m2, m_p1, m_p2
      DOUBLE PRECISION x_m1, x_m2, x_p1, x_p2
      DOUBLE PRECISION y_m1, y_m2, y_p1, y_p2
      DOUBLE PRECISION dxv
      INTEGER i

**************************************************************************
*  Definitions:
*      dx, dy = x(i+1) - x(i), y(i+1) - y(i)
*      m = slopes
**************************************************************************

**************************************************************************
*   Initializations
**************************************************************************

      DO 10 i = 1, aspline_size
         dx(i) = 0.
         dy(i) = 0.
         m(i) = 0.

         b(i) = 0.
         c(i) = 0.
         d(i) = 0.
10    CONTINUE

************************************************************************
*  Calculate spline values and store in "xfit", "yfit"
*    o  first, we make sure that there is data (ie., n_data > 0)
************************************************************************

      IF (n_data .GT. 0) THEN

************************************************************************
*  Here is the Akima spline fitting algorithm
************************************************************************

* get the endpoint slope calculations out of the way

         dx(1) = x(2) - x(1)
         dy(1) = y(2) - y(1)
         m(1) = dy(1) / dx(1)

         dx(2) = x(3) - x(2)   
         dy(2) = y(3) - y(2)
         m(2) = dy(2) / dx(2)

         dx(3) = x(4) - x(3)   
         dy(3) = y(4) - y(3)
         m(3) = dy(3) / dx(3)

         dx(n_data-1) = x(n_data) - x(n_data-1)
         dy(n_data-1) = y(n_data) - y(n_data-1)
         m(n_data-1) = dy(n_data-1) / dx(n_data-1)
             
         dx(n_data-2) = x(n_data-1) - x(n_data-2)
         dy(n_data-2) = y(n_data-1) - y(n_data-2)
         m(n_data-2) = dy(n_data-2) / dx(n_data-2)

         dx(n_data-3) = x(n_data-2) - x(n_data-3)
         dy(n_data-3) = y(n_data-2) - y(n_data-3)
         m(n_data-3) = dy(n_data-3) / dx(n_data-3)

* interpolate missing points

         x_m1 = x(1) + x(2) - x(3)
         y_m1 = (x(1)-x_m1) * (m(2) - 2. * m(1)) + y(1)

         m_m1 = (y(1)-y_m1)/(x(1)-x_m1)
   
         x_m2 = 2.d0 * x(1) - x(3)
         y_m2 = (x_m1-x_m2) * (m(1) - 2. * m_m1) + y_m1
   
         m_m2 = (y_m1-y_m2) / (x_m1-x_m2)

         x_p1 = x(n_data) + x(n_data-1) - x(n_data-2)
         y_p1 = (2.d0 * m(n_data-1) - m(n_data-2)) 
     +         * (x_p1 - x(n_data)) + y(n_data)

         m_p1 = (y_p1 - y(n_data)) / (x_p1 - x(n_data))
  
         x_p2 = 2.d0 * x(n_data) - x(n_data-2)
         y_p2 = (2.d0 * m_p1 - m(n_data-1)) * (x_p2 - x_p1) + y_p1
      
         m_p2 = (y_p2 - y_p1) / (x_p2 - x_p1)

* index = 1
         num = dabs(m(2) - m(1)) * m_m1 + dabs(m_m1 - m_m2) * m(1)
         den = dabs(m(2) - m(1)) + dabs(m_m1 - m_m2)
    
         IF (den .NE. 0.d0) THEN
            b(1) = num / den
         ELSE
            b(1) = 0.d0
         ENDIF

* index = 2 
         num = dabs(m(3) - m(2)) * m(1) + dabs(m(1) - m_m1) * m(2)
         den = dabs(m(3) - m(2)) + dabs(m(1) - m_m1)

         IF (den .NE. 0.d0) THEN
             b(2) = num / den
         ELSE
             b(2) = 0.d0
         ENDIF

* index = 3 -> n_data-1

         i = 0
         DO 30 i = 3, n_data-1

            dx(i) = x(i+1) - x(i)
            dy(i) = y(i+1) - y(i)
            m(i) = dy(i) / dx(i)

            dx(i+1) = x(i+2) - x(i+1)
            dy(i+1) = y(i+2) - y(i+1)
            m(i+1) = dy(i+1) / dx(i+1)


            num = dabs(m(i+1) - m(i)) * m(i-1) + 
     +            dabs(m(i-1) - m(i-2)) * m(i)
            den = dabs(m(i+1) - m(i)) + dabs(m(i-1) - m(i-2))

            IF (den .NE. 0.d0) THEN
                b(i) = num / den
            ELSE
                b(i) = 0.d0
            ENDIF
30       CONTINUE

* index = n_data - 1

         num = dabs(m_p1 - m(n_data-1)) * m(n_data-2)
     +        + dabs(m(n_data-2) - m(n_data-3)) * m(n_data-1)
         den = dabs(m_p1 - m(n_data-1)) 
     +        + dabs(m(n_data-2) - m(n_data-3))

         IF (den .NE. 0.d0) THEN
            b(n_data-1) = num / den
         ELSE
            b(n_data-1) = 0.d0
         ENDIF
 
* index = n_data

         num = dabs(m_p2 - m_p1) * m(n_data-1) +
     +         dabs(m(n_data-1) - m(n_data-2)) * m_p1
         den = dabs(m_p2 - m_p1) + dabs(m(n_data-1) - m(n_data-2))

         IF (den .NE. 0.d0) THEN
            b(n_data) = num / den
         ELSE 
            b(n_data) = 0.d0
         ENDIF

* assign Akima spline constants

         DO 40 i = 1, n_data
            
            dx(i) = x(i+1) - x(i)
            dy(i) = y(i+1) - y(i)
            m(i) = dy(i) / dx(i)
                      
            dxv = dx(i)
            c(i) = (3.d0 * m(i) - 2.d0 * b(i) - b(i+1)) / dxv
            d(i) = (b(i) + b(i+1) - 2.d0 * m(i)) / (dxv * dxv)

40       CONTINUE
      ENDIF          ! if n > 0

      RETURN
      END
