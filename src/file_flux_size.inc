*************************************************************************
*  Filename: file_flux_size.inc
*************************************************************************
*
*  Contains flux size for both 'sn_flux' arrays for models where data
*  is provided in files (replaces 'flux_size.inc').
*
*************************************************************************

      INTEGER file_time_size, file_energy_size, file_flux_size

***************************************************************************
* Bruenn 15 Ms (2000)
***************************************************************************

c      PARAMETER (file_time_size = 440)
c      PARAMETER (file_energy_size = 20)
c      PARAMETER (file_flux_size = file_time_size * file_energy_size)

***************************************************************************
* Bruenn 25 Ms --> black hole (2000)
***************************************************************************

c      PARAMETER (file_time_size = 190)
c      PARAMETER (file_energy_size = 20)
c      PARAMETER (file_flux_size = file_time_size * file_energy_size)

***************************************************************************
* Mezzacappa 13 Ms (2000)
***************************************************************************

      PARAMETER (file_time_size = 665)
      PARAMETER (file_energy_size = 12)
      PARAMETER (file_flux_size = file_time_size * file_energy_size)

