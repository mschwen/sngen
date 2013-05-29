***********************************************************************
*  Filename: sn_output.f
***********************************************************************
*
*  Generates SNOMAN output files with the SN neutrino time, energy
*  and class of event interaction 
*
*  Inputs:
*	o  array of interaction id number(s)
*	o  array of final state particle count totals for all specified
*	   interaction(s)
*
*************************************************************************

      SUBROUTINE sn_output(interaction_id,icount_part)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'

C  Common Block Declarations

      COMMON           /mcpl_stuff/
     +                 mcpl_flag, mcpl_file, mcpl_default
      INTEGER          mcpl_flag
      CHARACTER*40     mcpl_file, mcpl_default

      COMMON           /cmd_stuff/
     +                 cmd_flag, cmd_default, 
     +                 zdab_flag, zdab_default, hbook_default
      INTEGER          cmd_flag, zdab_flag
      CHARACTER*40     cmd_default, zdab_default, hbook_default

C  Global Variable Declarations

      INTEGER          interaction_id(*),   ! input
     +                 icount_part(*)      ! input

C  Local Variable Declarations

      INTEGER          i, interaction, num_counts
      LOGICAL          done

************************************************************************

************************************************************************
*  Check to see if we want to produce SNOMAN output of any sort
************************************************************************

      IF (mcpl_flag .EQ. 0 .AND. cmd_flag .EQ. 0) THEN
         WRITE (*,*)
         WRITE (*,*) 'Skipping SNOMAN output...'

         RETURN
      ELSE

************************************************************************
*  Initializations
************************************************************************

      i = 1
      done = .false.

************************************************************************

      WRITE (*,*)
      WRITE (*,*) 'Writing SNOMAN output:'

      DO WHILE (.NOT. done)

         interaction = interaction_id(i)

         IF (interaction .NE. 0) THEN

            num_counts = icount_part(interaction)

            IF (num_counts .GT. 0) THEN

************************************************************************
*  Output SNOMAN MCPL file
************************************************************************

               IF (mcpl_flag .EQ. 1) THEN
                  WRITE (*,*) 'Writing final particle MCPL files...'

                  CALL create_mcpl(interaction,num_counts)

               ENDIF

************************************************************************
*  Output SNOMAN .cmd file
************************************************************************

               IF (cmd_flag .EQ. 1) THEN
                  WRITE (*,*) 'Writing final particle cmd files...'

                  CALL create_cmd(interaction,num_counts)

               ENDIF

************************************************************************

            ENDIF    ! number of counts > 0

            i = i + 1

            done = .false.
         ELSE
            done = .true.
         ENDIF
      END DO

************************************************************************

      WRITE (*,*) 'All the SNOMAN output files have been created.'

      ENDIF

      RETURN
      END
