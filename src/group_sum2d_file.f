*************************************************************************
*  Filename: group_sum2d_file.f
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
 
      SUBROUTINE group_sum2d_file(interaction_nu,count_stat,t,filename)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'time_size.inc'

C  Global Variable Declarations

      INTEGER interaction_nu(*)                 ! input
      REAL*8 count_stat(interaction_total,*)    ! input
      REAL*8 t(*)                               ! input
      CHARACTER*40 filename                     ! input
      CHARACTER*40 file_dum

C  Common Block Declarations

      COMMON           /sn_time_entries/
     +                 time_entries1, time_entries2, time_entries
      INTEGER          time_entries1, time_entries2, time_entries 

C  Local Variable Declarations

      INTEGER interaction

      REAL*8 sum_count_stat(time_d_size)

      INTEGER i, j, k

*********************************************************************

*********************************************************************
*  Initializations
*********************************************************************

      DO i = 1, time_d_size
         sum_count_stat(i) = 0.d0
      END DO

      file_dum = filename
      filename = 'stats/'//file_dum

      OPEN (unit = 10, file = filename, status = 'unknown')

      file_dum = ''
      filename = ''
*********************************************************************
*  Read in count stats from the various interaction arrays
*********************************************************************

      i = 1

c      WRITE (*,*)
c      WRITE (*,*) 'Summing groups of interactions...'

      DO WHILE (      interaction_nu(i) .GT. 0 
     +          .AND. i .LE. interaction_total
     +          .AND. interaction_nu(i) .LE. interaction_total)

         interaction = interaction_nu(i)

         DO 10 j = 1, time_d_size
            sum_count_stat(j) = sum_count_stat(j) 
     +                         + count_stat(interaction,j)
10       CONTINUE

         i = i + 1
      END DO

* do not write zeros that occur for undelayed time
      DO 20 k = 1, time_d_size
         IF (k .LE. time_entries) WRITE (10,*) t(k), sum_count_stat(k)
         IF (k .GT. time_entries .AND. t(k) .NE. 0.)
     +    WRITE (10,*) t(k), sum_count_stat(k)
20    CONTINUE

      CLOSE (unit = 10)

      RETURN
      END
