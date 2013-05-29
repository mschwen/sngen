***************************************************************************
*  Filename = count_file_close.f
***************************************************************************

      SUBROUTINE count_file_close()

      IMPLICIT none

      COMMON           /sn_lumen_fit_file/
     +                 lumen_fit_file_flag,
     +                 lumen_nu_e_file,
     +                 lumen_nu_eb_file,
     +                 lumen_nu_mu_file
      INTEGER          lumen_fit_file_flag
      CHARACTER*40     lumen_nu_e_file,
     +                 lumen_nu_eb_file,
     +                 lumen_nu_mu_file

      COMMON           /FD_shape/
     +                 FD_time, FD_shape_file_flag,
     +                 FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu
      INTEGER          FD_shape_file_flag
      REAL*8           FD_time
      CHARACTER*40     FD_file_nu_e, FD_file_nu_eb, FD_file_nu_mu

      COMMON           /central_stuff/
     +                 central_flag, central_file
      INTEGER          central_flag
      CHARACTER*40     central_file

**************************************************************************
*  Check to see if we need to close files for the luminosity and energy
*  distributions
**************************************************************************

      IF (lumen_fit_file_flag .EQ. 1) THEN
         CLOSE (unit = 60)
         CLOSE (unit = 61)
         CLOSE (unit = 62)
      ENDIF

**************************************************************************
*  Check to see if we need to close files for the Fermi-Dirac distribution
**************************************************************************

      IF (FD_shape_file_flag .EQ. 1) THEN
         CLOSE (unit = 70)
         CLOSE (unit = 71)
         CLOSE (unit = 72)
      ENDIF

**************************************************************************
*  Check to see if we need to open files for saving the central values
**************************************************************************

      IF (central_flag .EQ. 1) THEN
         CLOSE (unit = 75)
      ENDIF

      RETURN
      END
