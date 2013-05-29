*************************************************************************
*  Filename: group_sum2d.f
*************************************************************************
*
*  This subroutine will find the total count statistic (countrate, 
*  cumulative counts, counts) for the specified set of interactions 
*  (for example: neutrino type, interaction type, region type).
*          
*  Inputs: 
*    	o  array of interaction id numbers for a particular class of
*          events 
*  Outputs:
*	o  array of count totals for all specified interaction(s)
*
*************************************************************************
 
      SUBROUTINE group_sum2d(interaction_nu,count_stat,sum_count_stat)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'sample_count_size.inc'
      INCLUDE 'sn_classes.inc'

C  Global Variable Declarations

      INTEGER interaction_nu(*)                 ! input
      REAL*8 count_stat(interaction_total,*)    ! input
      REAL*8 sum_count_stat                     ! output

C  Local Variable Declarations

      INTEGER i, j, interaction

*********************************************************************

*********************************************************************
*  Initializations
*********************************************************************

      sum_count_stat = 0.d0

*********************************************************************
*  Read in count stats from the various interaction arrays
*********************************************************************

      i = 1

c      WRITE (*,*)
c      WRITE (*,*) 'Summing groups of interactions...'

      DO WHILE (      interaction_nu(i) .GT. 0 
     +          .AND. i .LE. interaction_total + 1
     +          .AND. interaction_nu(i) .LE. last_HALO + 1)

         interaction = interaction_nu(i)

         DO 10 j = 1, sample_count_size
            sum_count_stat = sum_count_stat 
     +                         + count_stat(interaction,j)
10       CONTINUE
         i = i + 1
      END DO

      RETURN
      END
