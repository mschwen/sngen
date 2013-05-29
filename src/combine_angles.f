***************************************************************************
*  Filename: combine_angles.f
***************************************************************************
*
*  Given an initial set of directional cosines, U,V,W, a polar angle with
*  respect to those cosines, and the azimuthal angle, this routine will
*  calculate the new directional cosines, and return them in the 
*  array: 'part_dir_cos'
*
*   Inputs:
*	o  supernova directional cosines
*	o  scattering angle with respect to supernova directional cosines
*	o  azimuthal angle
*   Outputs:
*	o  sampled final state particle directional cosines 
*          (combined angle directional cosines)
*
***************************************************************************

      SUBROUTINE combine_angles(sn_dir_cos,scat_angle,azimuth_angle, 
     +                          part_dir_cos)
   
      IMPLICIT none

C   Global Variable Declarations

      REAL*8 sn_dir_cos(*)
      REAL*8 scat_angle, azimuth_angle
      REAL*8 part_dir_cos(3)

C   Local Variable Declarations

      REAL*8 u, v, w

      REAL*8 renorm
      REAL*8 c_angle, s_angle
      REAL*8 c_azimuth, s_azimuth
      REAL*8 s_theta
      REAL*8 c_phi, s_phi

***************************************************************************

***************************************************************************
*  Initializations
***************************************************************************

**** set supernova direction cosines calculated earlier

      u = sn_dir_cos(1)
      v = sn_dir_cos(2)
      w = sn_dir_cos(3)

****************************************************************************
*  Assign new angles based on scattering about the supernova initial
*  direction (taken from the SNOMAN routine 'rayscatter.for')
****************************************************************************

***  Make sure directional cosines are normalized

      renorm = dsqrt(u*u + v*v + w*w)
      u = u / renorm
      v = v / renorm
      w = w / renorm

***  Define useful trig functions

      c_angle = dcos(scat_angle)
      s_angle = dsin(scat_angle) 

      c_azimuth = dcos(azimuth_angle)
      s_azimuth = dsin(azimuth_angle) 
 
*****************************************************************************
*  For an explanation of the geometrical calculations, see the SNOMAN
*  manual.  The angle will be between 0 and 180 degrees, so the sine
*  will always be positive hence the 'sqrt()' will give the correct answer
******************************************************************************

      IF (w .GE. 1.d0) THEN
        s_theta = 0.d0
      ELSE
        s_theta = dsqrt(1.d0 - w*w)
      ENDIF

***  If sin(theta)=0, then phi is not defined. This must be error trapped.

      IF (s_theta .GT. 0.00001)THEN
        c_phi = u / s_theta
        s_phi = v / s_theta
      ELSE
        c_phi = 1.d0
        s_phi = 0.d0
      ENDIF
 
***  Calculate the particle's new combined angle directional cosines

      part_dir_cos(1) = c_angle*u + s_angle*(c_azimuth*w*c_phi 
     +                  - s_azimuth*s_phi)
      part_dir_cos(2) = c_angle*v + s_angle*(c_azimuth*w*s_phi 
     +                  + s_azimuth*c_phi)
      part_dir_cos(3) = c_angle*w - s_angle*c_azimuth*s_theta

      IF (     part_dir_cos(1) .GT. 1.d0
     +    .OR. part_dir_cos(1) .LT. -1.d0) THEN
         WRITE (*,*) '*** combine_angles.f: Something went wrong!'
         WRITE (*,*) '*** combine_angles.f: dir_cos(1) =',
     +                part_dir_cos(1)
      ENDIF

***************************************************************************

      RETURN
      END
