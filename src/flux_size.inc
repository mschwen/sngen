*************************************************************************
*  Filename: flux_size.inc
*************************************************************************
*
*  Contains flux size for both 'sn_flux' and 'det_flux' arrays
*
*    o time_size   = (tstop - tstart) / tstep + 1
*    o energy_size = (emax - emin) / estep
*
*    -> flux_size = time_size * energy_size
*
*************************************************************************

      INTEGER flux_size

*************************************************************************
*  Burrows et al. (1992)
*************************************************************************

* time: -100ms ->  1000ms; 1ms   energy: 0.1 -> 100MeV; 0.1MeV = 1101000
c* time: 1010ms -> 50000ms; 10ms  energy: 0.1 -> 100MeV; 0.1MeV = 4900000
* time: 1020ms -> 50000ms; 20ms  energy: 0.1 -> 100MeV; 0.1MeV = 2450000

c      PARAMETER (flux_size = 6001000)   ! Burrows default
c      PARAMETER (flux_size = 3551000)   ! Burrows default (sngen default)

*************************************************************************
*  Beacom & Vogel (1998)
*************************************************************************

c* time: -90ms -> 20000ms; 5ms   energy: 0.1 -> 170MeV; 0.1MeV = 6832300
* time: -90ms -> 20000ms; 10ms   energy: 0.1 -> 170MeV; 0.1MeV = 3417000

c      PARAMETER (flux_size = 6832300)  ! BeacomVogel default
      PARAMETER (flux_size = 3417000)  ! BeacomVogel default (< Burrows)

***************************************************************************
* Bruenn, Mezzacappa (2001) --> see file_flux_size
***************************************************************************
