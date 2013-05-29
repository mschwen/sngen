**************************************************************************
*  Filename: poisson.f
**************************************************************************
*
*  Deals with Poisson distribution of detection process
*
*  Inputs:
*       o  this interaction id (not array)
*	o  average number of detector counts from "count_sum" routine
*  Outputs:
* 	o  whole number of detector counts distributed according
*          to a Poisson distribution
*
**************************************************************************

      SUBROUTINE poisson(this_interaction,count_ave,count_pois)

      IMPLICIT none

      INCLUDE 'random_seed.inc'

C  Common Block Declarations

      COMMON           /sn_pois_flag/
     +                 pois_flag
      INTEGER          pois_flag

C  Global Variable Declarations

      INTEGER          this_interaction    ! input
      REAL*8           count_ave            ! input
      REAL*8           count_pois           ! output

C  Local Variable Declarations

      REAL*8           poidev  ! function

**************************************************************************

      IF (pois_flag .EQ. 1) THEN

         WRITE (*,*) '  Applying Poisson statistics...'

         IF (count_ave .GT. 0.) THEN

            count_pois = poidev(count_ave,iseed)

            WRITE (*,5) count_ave, count_pois
5           FORMAT (3x,'-->', f7.2,' counts becomes ',f7.2,' counts')
         ENDIF

      ELSEIF (pois_flag .EQ. 0) THEN

         WRITE (*,*) '  Skipping Poisson statistics...' 

         count_pois = count_ave
      ENDIF

      RETURN
      END
