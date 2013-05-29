***************************************************************************
*  Filename: xsect_nu_d_diff.f
***************************************************************************
* 
* o   Just re-maps input angle and energy into a function for use in
*     minimization (Powell's method).
*
*     Same as SNOMAN routine: "f_diff_ns.for" for NSAPMGK cross sections
*
*     Program Notes:-
*     =============
* 
*     If cosine of selected angle of outgoing particle is outside physical
*     boundaries (+/-1) or selected energy is negative, then returns a large
*     positive value to ensure that minimum will be found in physical region.
*
*****************************************************************************
 
      REAL FUNCTION xsect_nu_d_diff(p)

      IMPLICIT none

      INCLUDE 'physics_constants.inc'
 
C   Common Block Declarations

      REAL E_nu
      INTEGER int_type
      COMMON /nu_d_en_stuff/ E_nu, int_type ! need to do this for 'powell'

C   Global Variable Declarations

      REAL p(*)

C   Local Variable Declarations

      REAL cos_theta, E_part
      REAL f_xs_diff

*************************************************************************** 

      cos_theta = p(1)
      E_part = p(2) ! KN cross sections require kinetic E
c      E_part = p(2) + m_electron  ! new NSAPMGK cross sections require total E

      IF(     (abs(cos_theta) .GT. 1.d0)
     +   .OR. (E_part .LT. 0.d0)        ) THEN
	xsect_nu_d_diff = 1.d6      !force physical boundaries
	RETURN
      ELSE
         xsect_nu_d_diff = -f_xs_diff(int_type,E_nu,cos_theta,E_part)
      ENDIF

      RETURN
      END
