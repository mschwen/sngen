*************************************************************************
*  Filename: group_isum1d.f
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
 
      SUBROUTINE group_isum1d(interaction_nu,icounts,sum_icounts)

      IMPLICIT none

      INCLUDE 'sn_classes.inc'

C  Global Variable Declarations

      INTEGER          interaction_nu(*),     ! input
     +                 icounts(*),            ! input
     +                 sum_icounts            ! output

C  Local Variable Declarations

      INTEGER          i, interaction

*********************************************************************

*********************************************************************
*  Initializations
*********************************************************************

      sum_icounts = 0

*********************************************************************
*  Read in count stats from the various interaction arrays
*********************************************************************

      i = 1

c      WRITE (*,*)
c      WRITE (*,*) 'Summing groups of interactions...'

      DO WHILE (      interaction_nu(i) .GT. 0 
     +          .AND. i .LE. last_HALO + 1
     +          .AND. interaction_nu(i) .LE. last_HALO + 1)

         interaction = interaction_nu(i)
         sum_icounts = sum_icounts + icounts(interaction)

         i = i + 1
      END DO

      RETURN
      END
