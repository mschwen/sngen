************************************************************************
*  Filename: sngen.f
************************************************************************
*
*  Initial outline: Siong Ng
*  Complete Rework/Revisions: Jaret Heise
*
*  This is the main body of the SNO supernova generator program. 
*  All of the user inputs are to be found in the "sn_user.dat"
*  file, including distances, neutrino masses, supernova direction, 
*  output data files and all flags.
*
*************************************************************************

      PROGRAM sngen

      IMPLICIT none

      INCLUDE 'version.inc'
      INCLUDE 'interaction_total.inc'

C  Local Variable Declarations

      INTEGER          i
      INTEGER          id_interaction(100), 
     +                 this_interaction,
     +                 exp_id


      REAL*8           this_count_ave,
     +                 this_count_pois, count_pois(100),
     +                 this_count_part, count_part(100),
     +                 sum_neutrino_count, Test
      INTEGER          icount_pois(100),
     +                 this_icount_part, icount_part(100),
     +                 n_particles, id_part(3)

      DATA             sum_neutrino_count  /0.d0/,
     +                 count_pois          /100 * 0.d0/,
     +                 count_part          /100 * 0.d0/,
     +                 icount_pois         /100 * 0/,
     +                 icount_part         /100 * 0/

************************************************************************

      WRITE (*,5) major_version, minor_version
5     FORMAT ('sngen version ',i1,'.',i1)
      WRITE (*,*)

************************************************************************
*  Read SN user inputs into COMMON blocks from file: "dat/sn_user.dat"
************************************************************************

      CALL read_user()

************************************************************************
* Choose experiment
************************************************************************

       WRITE(*,*) "Choose experiment (SNO = 1, SNOPLUS = 2, HALO = 3):"
       READ(*,*) exp_id

************************************************************************
*  Select supernova neutrino reaction(s) to study
************************************************************************

      CALL sn_interactions(id_interaction,exp_id)

************************************************************************
* Write neutrino flux to file (specifically to verify flavour swapping)
* This program was also used to calculate the critical energies during
* flavour swapping 
* If wanting to do this, the user should only choose one interaction,
* this call should be commented for normal operation
************************************************************************

c      CALL flux_integral(id_interaction)

************************************************************************
*  Loop over specified supernova interactions
************************************************************************

      CALL count_file_open()

      OPEN (unit=9, file='sum_total.out', status='unknown')

      DO 10 i = 1, interaction_total

* analyze specified interactions

         IF (id_interaction(i) .GT. 0) THEN
            this_interaction = id_interaction(i)

************************************************************************
*  Integrate supernova source and detector fluxes to get both the
*  supernova neutrino count rate as well as the total number of counts
*  for the specified interactions
************************************************************************

            CALL count_sum(this_interaction,this_count_ave)

            sum_neutrino_count = sum_neutrino_count + this_count_ave

***********************************************************************
*  Distribute detector counts according to Poisson distribution
***********************************************************************

            CALL poisson(this_interaction,this_count_ave,
     +                   this_count_pois)

            WRITE (9,*) this_interaction, this_count_pois

            count_pois(this_interaction) = this_count_pois
            icount_pois(this_interaction) = nint(this_count_pois)
         
************************************************************************
*  Sample from the batch of neutrinos that were created and produce
*  the number indicated by the count_rate rountine and then sample the
*  final state products that SNO could possibly detect
************************************************************************

            CALL sample(this_interaction,this_count_pois,
     +                  interaction_total)

            CALL get_num_part(this_interaction,n_particles,id_part)
            this_count_part = this_count_pois * dble(n_particles)
            this_icount_part = nint(this_count_pois) * n_particles

            count_part(this_interaction) = this_count_part
            icount_part(this_interaction) = this_icount_part
         ENDIF
10    CONTINUE

      CALL count_file_close()

      CLOSE (unit = 9)

      WRITE (*,*)
      WRITE (*,*) 'All the specified interactions have been generated an
     +d sampled'

      WRITE (*,*)
      WRITE (*,*) 'The total average number of SN neutrinos is',
     +             sum_neutrino_count

*************************************************************************
*  Produce a summary file
*************************************************************************
c      CALL flux_integral(id_interaction)

      CALL summary(id_interaction,count_pois,icount_pois,
     +             count_part,icount_part,exp_id)

*************************************************************************
*  Produce some additional statistics
*************************************************************************

      CALL stats(exp_id,id_interaction,count_pois,icount_pois)

************************************************************************
*  Simulate SNO data acquisition trigger gate in order to see if there
*  are any events that overlap in time (currently, just checks neutrino times)
************************************************************************

      CALL sn_daq(id_interaction,count_pois)

************************************************************************
*  Output data so that SNOMAN can read it
************************************************************************

      CALL sn_output(id_interaction,icount_part)

************************************************************************

      WRITE (*,*)
      WRITE (*,*) 'C''est Fini!'
	
      STOP
      END
