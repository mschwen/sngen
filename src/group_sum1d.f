*************************************************************************
*  Filename: group_sum1d.f
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
 
      SUBROUTINE group_sum1d(interaction_nu,counts,sum_counts,exp_id)

      IMPLICIT none

      INCLUDE 'sn_classes.inc'
      INCLUDE 'interaction_total.inc'

C  Global Variable Declarations

      INTEGER          interaction_nu(*)     ! input
      INTEGER          exp_id                ! input
      REAL*8           counts(*)             ! input
      REAL*8           sum_counts            ! output

C  Local Variable Declarations

      INTEGER i, interaction

*********************************************************************

*********************************************************************
*  Initializations
*********************************************************************

      sum_counts = 0.d0

*********************************************************************
*  Read in count stats from the various interaction arrays
*********************************************************************

      i = 1

c      WRITE (*,*)
c      WRITE (*,*) 'Summing groups of interactions...'

         DO WHILE (      interaction_nu(i) .GT. 0
     +          .AND. i .LT. interaction_total + 1 
     +          .AND. interaction_nu(i) .LT. last_HALO + 1)

         interaction = interaction_nu(i)
         sum_counts = sum_counts + counts(interaction)

c         WRITE (*,*) 'i, interaction_nu(i), interaction, counts = ',
c     +                i, interaction_nu(i), interaction, 
c     +                counts(interaction)
         i = i + 1
         END DO

      RETURN
      END
