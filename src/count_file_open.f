*****************************************************************************
*  Filename: count_file_open.f
*****************************************************************************

      SUBROUTINE count_file_open()

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
*  Check to see if we need to open files for the luminosity and energy
*  distributions (open them here, but write to them in the 'sn_num'
*  function)
**************************************************************************

      IF (lumen_fit_file_flag .EQ. 1) THEN
         OPEN (unit=60,file=lumen_nu_e_file,status='unknown')
         OPEN (unit=61,file=lumen_nu_eb_file,status='unknown')
         OPEN (unit=62,file=lumen_nu_mu_file,status='unknown')
      ENDIF

**************************************************************************
*  Check to see if we need to open files for the Fermi-Dirac distribution
*  (open them here, but write to them in the 'fd_spectrum' function)
**************************************************************************

      IF (FD_shape_file_flag .EQ. 1) THEN
         OPEN (unit=70,file=FD_file_nu_e,status='unknown')
         OPEN (unit=71,file=FD_file_nu_eb,status='unknown')
         OPEN (unit=72,file=FD_file_nu_mu,status='unknown')
      ENDIF

**************************************************************************
*  Check to see if we need to open files for saving the central values
**************************************************************************

      IF (central_flag .EQ. 1) THEN
         OPEN (unit=75,file=central_file,status='unknown')
      ENDIF

      RETURN
      END
