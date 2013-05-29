**************************************************************************
*  Filename: det_flux.f
**************************************************************************
*
*  This routine will return the supernova flux expected in the SNO
*  detector at a given energy and time for a given supernova interaction.
*  Time is not explicitly required in this routine.
*
*  Calculations performed in this routine (or below) include:
*	o  interaction cross-sections with SNO detector volumes
*
*  Inputs:
*       o  supernova interaction id number  (not array) 
*	o  supernova neutrino energy  [MeV]
*	o  supernova source flux at specified distance  [#/(s*MeV*cm^2)]
*  Outputs:
*	o  supernova flux in the SNO detector  [#/(s*MeV)]
*
*  Note on units: Detector flux becomes [#/(s*MeV*m^2)] after scaling for
*                 SN distance but the cross-section dimensions [m^2] 
*                 cancel this additional factor later
*
**************************************************************************

      FUNCTION det_flux(interaction,E_nu,sn_flux)

      IMPLICIT none

      INCLUDE 'trig_constants.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON            /distance/
     +                  distance
      REAL*8            distance

C  Global Variable Declarations

      INTEGER           interaction           ! input
      REAL*8            E_nu, sn_flux         ! input

C  Local Variable Declarations

      REAL*8            det_flux     ! function
      REAL*8            normalize, N_target, xsect_total

**************************************************************************

* normalize supernova flux to the supernova distance to the Earth
      
      normalize = 4.d0 * pi * (distance * kpc2m)**2

      sn_flux = sn_flux / normalize

* get the number of target atoms from the appropriate interaction region

      CALL get_target(interaction,N_target)

* get cross section [m^2]

      CALL get_xsect(interaction,E_nu,xsect_total)

* compute the supernova neutrino detector flux [#/(MeV*s)]

c      WRITE (*,*) E_nu, xsect_total, sn_flux


      det_flux = N_target * xsect_total * sn_flux

      END
