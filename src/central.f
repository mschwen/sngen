***********************************************************************
*  Filename: central.f
***********************************************************************
*
*  Calculate central values as defined by Beacom.
*
***********************************************************************

      FUNCTION central(this_interaction,E_nu,time)
   
      IMPLICIT none

C   Global Variable Declarations

      INTEGER this_interaction           ! input
      REAL*8 E_nu, time                  ! input 

C   Local Variable Declarations

      REAL*8 central, fd_spectrum        ! functions
      REAL*8 xsect_total
      INTEGER flux_index

*************************************************************************

* call this routine to assign "this_energy_index"
* (passed through COMMON block)

      CALL get_flux_index(E_nu,time,flux_index)

* get total cross section

      CALL get_xsect(this_interaction,E_nu,xsect_total)

* calculate central values

      central = xsect_total * fd_spectrum(this_interaction,E_nu,time)

c      write (*,*) this_interaction,E_nu,time,xsect_total,
c     +            fd_spectrum(this_interaction,E_nu,time),
c     +            central

      END
