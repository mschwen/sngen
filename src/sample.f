***********************************************************************
*  Filename: sample.f
***********************************************************************
*
*  Samples supernova neutrino distribution to give time and energy of
*  Poisson distributed events
*
*  Inputs:
*	o  this interaction id (not array)
*	o  Poisson count total for this interaction
*
*  Outputs:
*	o  sampled neutrino counts (output file from "sample_nu")
*	o  sampled final state particle counts (output file from "sample_part")
*
*************************************************************************

      SUBROUTINE sample(this_interaction,count_pois)

      IMPLICIT none

C  Global Variable Declarations

      INTEGER          this_interaction            ! input
      REAL*8           count_pois                  ! input

C  Common Block Declarations

      COMMON           /generator_flag/
     +                 gen_flag, counts_file
      INTEGER          gen_flag
      CHARACTER*40     counts_file

      COMMON           /number_supernova/
     +                 n_supernova
      INTEGER          n_supernova

      COMMON           /nu_sample_file/
     +                 nu_sample_file
      CHARACTER*40     nu_sample_file

      COMMON           /part_sample_file/
     +                 part_sample_file
      CHARACTER*40     part_sample_file

C  Local Variable Declarations

      INTEGER          file_flag
      CHARACTER*40     sample_file
      INTEGER          i, j, icount_neutrino, icount_particle,
     +                 n_particles, id_part

************************************************************************

* check to see if we need to sample
* (sample files may already exist and all we need to do is generate the
*  SNOMAN files)

      IF (     gen_flag .EQ. 2
     +    .OR. gen_flag .EQ. 3) THEN
         WRITE (*,*) '  Skipping all sampling...'      
         RETURN
      ENDIF

* check to see if we want to load in count numbers
* (in order to get statistics with a certain number of counts on all reactions)

      IF (gen_flag .EQ. 1) THEN

         OPEN (unit=25, file=counts_file, status='unknown')

* read counts file until the line with the specified interaction
* (count_ave should then have the correct counts for this interaction)

c         WRITE (*,*) '  Reading supernova neutrino counts for interactio
c     +n ', this_interaction,'...'

         DO i = 1, this_interaction
            READ (25,*) j, count_pois    ! interaction id, counts
         END DO
         CLOSE (unit = 25)

         count_pois = count_pois * dble(n_supernova)

c         WRITE (*,5) this_interaction, count_pois
c5        FORMAT (3x,'For interaction ',i2,' the total number of counts i
c     +s now ' ,f9.3)
      ENDIF

************************************************************************
*  Neutrino sampling
************************************************************************

      IF (this_interaction .GT. 0) THEN

         icount_neutrino = nint(count_pois)  ! round events to integer

* Check to see whether this sample file already exists or not

         file_flag = 0

         CALL file_increment(nu_sample_file,this_interaction,
     +                       sample_file)
         CALL file_check(sample_file,file_flag)

* Sampling file already exists

         IF (file_flag .EQ. 1) THEN

            WRITE (*,*) '  Sampling interaction', this_interaction,'...' 
            WRITE (*,*) '  -> ', icount_neutrino,
     +                  ' neutrinos have already been generated.'
            WRITE (*,*) '  -> skipping to next section...'

* Sampling file does not exist

         ELSEIF (file_flag .EQ. 0) THEN

**********************************************************************
*  Sample the SN neutrino distribution over time and energy
**********************************************************************

            IF (icount_neutrino .GT. 0) THEN
               CALL sample_nu(this_interaction,icount_neutrinO)
            ELSE
               WRITE (*,*) '  --> no neutrinos to sample'
            ENDIF

         ELSE
            WRITE (*,*) '*** sn_sample.f: Something went wrong!'
            WRITE (*,*) '*** sn_sample.f: (sample_nu section)'
         ENDIF  ! check if neutrino sample file has already been generated

*************************************************************************
*  Particle sampling
*************************************************************************

* Check to see whether this interaction file already exists or not

         file_flag = 0

         CALL file_increment(part_sample_file,this_interaction,
     +                       sample_file)
         CALL file_check(sample_file,file_flag)

* Particle sample file exists already

         IF (file_flag .EQ. 1) THEN

            CALL get_num_part(this_interaction,n_particles,id_part)

            icount_neutrino = nint(count_pois) ! have to do this again for
                                               ! some reason *boggle*
            icount_particle = icount_neutrino * n_particles

            WRITE (*,*) '  Sampling interaction', this_interaction,'...' 
            WRITE (*,*) '  -> ', icount_particle,
     +                  ' particles have already been generated.'
            WRITE (*,*) '  -> skipping to next section...'

* Particle sample file does not exist

         ELSEIF (file_flag .EQ. 0) THEN

**********************************************************************
*  Sample the SN neutrino distribution for time and energy
**********************************************************************

            IF (icount_neutrino .GT. 0) THEN
               CALL sample_part(this_interaction,icount_neutrino)
            ELSE
               WRITE (*,*) '  --> no particles to sample'
            ENDIF

         ELSE
            WRITE (*,*) '*** sn_sample.f: Something went wrong!'
            WRITE (*,*) '*** sn_sample.f: (sample_part section)'
         ENDIF     ! check if particle sample file has already been generated
      ENDIF        ! interaction id > 0

      RETURN
      END
