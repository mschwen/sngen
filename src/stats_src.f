***********************************************************************
*  Filename: stats_src.f
***********************************************************************
*
*  Inputs:
*	o  neutrino type 
*	o  neutrino emission time (ms)
*  Outputs:
*       o total SN average neutrino energy for the specified neutrino type
*
*  Note: this routine only depends on time, so when integrating, we should
*  save the previous number 
* 
***********************************************************************
     
      SUBROUTINE stats_src(nu_type,total_E,total_num,total_Eavg)

      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON           /time_limits/ 
     +                 time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2,
     +                 dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2
      REAL*8           time_start1, time_stop1, tstep1,
     +                 time_start2, time_stop2, tstep2
      INTEGER          dum_tmin1, dum_tmax1, dum_tstep1,
     +                 dum_tmin2, dum_tmax2, dum_tstep2

C  Global Variable Declarations

      INTEGER          nu_type                           ! input
      REAL*8           total_E, total_num, total_Eavg    ! output

C  Local Variable Declarations

      REAL*8           sn_lum, sn_Eavg  ! functions

      INTEGER          t_dum
      REAL*8           time

      REAL*8           t_fac ! factor to improve integration resolution
      PARAMETER        (t_fac = 1000.d0) ! gives 0.001 ms steps
      REAL*8           d_lum, d_Eavg, t_step, sum_lum, sum_num

************************************************************************

      sum_lum = 0.d0
      sum_num = 0.d0

      DO 10 t_dum = nint(time_start1*t_fac), nint(time_stop2*t_fac)

         time = float(t_dum) / t_fac  ! 1 ms steps (all times are in ms)
         t_step = 1.d0 / t_fac 

         d_lum  = sn_lum(nu_type,time)
         d_Eavg = sn_Eavg(nu_type,time)

* include contributions to the total energy from both muon and tau types
* (in the "stats.f" routine, nu_mu is taken to mean both muon and tau)
         IF (     nu_type .EQ. id_nu_mu
     +       .OR. nu_type .EQ. id_nu_mubar) THEN
            d_lum = d_lum * 2.d0
         ENDIF

         sum_lum = sum_lum + d_lum * t_step
         sum_num = sum_num + (d_lum / d_Eavg) * t_step

10    CONTINUE

      total_E    = sum_lum * ms2s
      total_num  = sum_num * ms2s
      total_Eavg = sum_lum / sum_num  ! unit conversions cancel out here

      END
