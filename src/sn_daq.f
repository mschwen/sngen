***********************************************************************
*  Filename: sn_daq.f
***********************************************************************
*
*  Simulation of the SNO daq trigger gate.  Basically, this routine
*  checks all neutrino events time to see if there is any overlap
*  on the order of the Global Trigger lockout time.
*
*  Currently, this only checks the times for neutrinos (as opposed to
*  final state particles).
*
*  Inputs:
*	o  array of interaction id numbers
*	o  array of Poisson count totals for all specified interactions
*
*************************************************************************

      SUBROUTINE sn_daq(interaction_id,count_pois)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'sn_count_size.inc'

C  Common Block Declarations

      COMMON /nu_sample_file/ nu_sample_file
      CHARACTER*40 nu_sample_file

      COMMON           /sn_time_overlap/ 
     +                 trigger_gate, 
     +                 time_overlap_flag, time_overlap_file
      REAL*8 trigger_gate
      INTEGER time_overlap_flag
      CHARACTER*40 time_overlap_file

C  Global Variable Declarations

      INTEGER interaction_id(*)
      REAL*8 count_pois(*)

C  Local Variable Declarations

      INTEGER kth_interaction, interaction
      INTEGER nevents
      INTEGER i, j, l

      REAL*8 time_d(sn_count_size), en(sn_count_size)

      INTEGER int_i, int_j
      INTEGER event_i, event_j

      REAL*8 nu_time(interaction_total,sn_count_size)
      INTEGER count, event_count(sn_count_size)
      REAL*8 time_diff

      INTEGER lost_trigger

************************************************************************

      IF (time_overlap_flag .EQ. 1) THEN

         WRITE (*,*)
         WRITE (*,*) 'Simulating electronics gate with time gate:',
     +                trigger_gate,'ms ...'

************************************************************************
*  Initializations
************************************************************************

         lost_trigger = 0

         DO 10 i = 1, interaction_total
            DO 11 j = 1, sn_count_size 
               nu_time(i,j) = -1.d4
11          CONTINUE
10       CONTINUE

***********************************************************************

         kth_interaction = 1

         DO WHILE (kth_interaction .LE. interaction_total)

            interaction = interaction_id(kth_interaction)

***  Round all interaction count values to integers

               nevents = nint(count_pois(interaction))

***********************************************************************
*  Read in all nu_sample files and build time array for comparison
***********************************************************************

               IF (nevents .GT. 0) THEN

                  CALL file_open(25,nu_sample_file,interaction)

***  Reset interaction counter for each interaction
                  count = 0

                  DO 21 l = 1, nevents
                     READ (25,*) time_d(l), en(l)

                     count = count + 1          ! per interaction
                     nu_time(interaction,count) = time_d(l)
21                CONTINUE

                  event_count(interaction) = count

                  CLOSE (unit = 25)
               ENDIF

               kth_interaction = kth_interaction + 1
         END DO

**********************************************************************
*  Output all event times to single file
**********************************************************************

c         OPEN (unit=18, file='event_list.out', status='unknown')

c         DO 25 int_i = 1, interaction_total
c            DO 26 event_i = 1, event_count(int_i)
c               WRITE (18,*) nu_time(int_i,event_i)
c26          CONTINUE
c25       CONTINUE

c         CLOSE (unit = 18)

**********************************************************************
*  Check for time overlap that will cause one of the neutrinos to be
*  'lost' as far as the global trigger is concerned
**********************************************************************

         OPEN (unit=19, file=time_overlap_file, status='unknown')

         DO 30 int_i = 1, interaction_total
           IF (event_count(int_i) .GT. 0) THEN
             DO 31 event_i = 1, event_count(int_i)
               DO 32 int_j = 1, interaction_total
***  Don't want to double count interaction comparisons
***  (eg. don't compare interaction 1 with interaction 2, then 2 with 1; 
***   interaction 1 with interaction 1 is fine)
                 IF (      (int_i .LE. int_j)
     +               .AND. (event_count(int_j) .GT. 0)) THEN
c                   WRITE (*,*) '*comparing int ',int_i,' with int',int_j
                   DO 33 event_j = 1, event_count(int_j)
***  Make sure we haven't gone past the number of events
***  (the array is set to 'sn_count_size' so there is room off the end)
                     IF (      (nu_time(int_i,event_i) .NE. -1.d4) 
     +                 .OR. (nu_time(int_j,event_j) .NE. -1.d4)) THEN
***  Don't want to double count event comparisons within an interaction
***  (eg. don't compare event 1 with event 2, then 2 with 1 in the same
***   interaction)
                       IF (      ((int_i .EQ. int_j)
     +                     .AND. (event_i .LT. event_j))
     +                      .OR. (int_i .NE. int_j)     ) THEN
c                         WRITE (*,*) 'comparing ev ', event_i,
c     +                    ' with ev ', event_j
c                        WRITE (*,*) nu_time(int_i,event_i)
c                        WRITE (*,*) nu_time(int_j,event_j)
                         time_diff = dabs(nu_time(int_i,event_i)
     +                              - nu_time(int_j,event_j))
                         IF (time_diff .LE. trigger_gate) THEN
                            lost_trigger = lost_trigger + 1
                            WRITE (19,*) int_i,event_i,int_j,event_j,
     +                                   time_diff

                         ENDIF  ! time_diff test
                       ENDIF    ! event double count test
                     ENDIF      ! array bounds test
33                 CONTINUE     ! jth_event loop
                 ENDIF          ! interaction double counting test
32             CONTINUE         ! jth_interaction loop
31           CONTINUE           ! ith_event loop
           ENDIF                ! non-zero counts test
30       CONTINUE               ! ith_interaction loop

         CLOSE (unit = 19)

         WRITE (*,*) 'The lost trigger count due to time overlap is',
     +                lost_trigger,'.'

      ELSE
         WRITE (*,*)
         WRITE (*,*) 'Skipping DAQ trigger section...'
      ENDIF

      RETURN
      END
