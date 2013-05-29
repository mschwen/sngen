***********************************************************************
*  Filename: sample_part.f
***********************************************************************
*
*  Given a neutrino time and energy, this routine will sample as many
*  final state particles as can be detected by the SNO detector.  The idea
*  is to generate a random energy for the detected particle between the
*  maximum and minimum kinematically allowed kinteic energies.  Knowing
*  the neutrino and detected particle energy the angular distribution can
*  then be evaluated.
*
*  Note: As input, SNOMAN takes _total energy_ for electrons/positrons
*        and _kinetic energy_ for neutrons. 
*
*  Inputs:
*	o  interaction number (not array)
*	o  number of sample events to produce
*	o  supernova direction cosines
*  Outputs:
*	o  number of final state particles for this interaction
*       o  output file ('part_sample_file')
*	     --> particle id, time, energy, direction(x,y,z)
*
*********************************************************************** 

      SUBROUTINE sample_part(interaction,n_events)

      IMPLICIT none

      INCLUDE 'sample_count_size.inc'
      INCLUDE 'id_int.inc'
      INCLUDE 'z.inc'
      INCLUDE 'bank_mctk.inc'
      INCLUDE 'titles_msw.inc'  ! for C_SCTN flag

C  Common Block Declarations

      COMMON            /part_sample_file/
     +                  part_sample_file
      CHARACTER*100      part_sample_file
      COMMON            /nu_sample_file/
     +                  nu_sample_file
      CHARACTER*100      nu_sample_file

C  Global Variable Declarations

      INTEGER           interaction,           ! input
     +                  n_events               ! input

C  Local Variable Declarations

      LOGICAL           first_call

      REAL*8            time_d(sample_count_size), 
     +                  E_nu(sample_count_size)
      REAL*8            decay_delay

      INTEGER           int_type, nu_type

      INTEGER           n_particles_total, n_particles, id_part(3)

      REAL*8            E_part
      REAL*8            sn_dir_cos(3), part_dir(3)

* nu_e_bar + d buffer variables
      INTEGER           redo, id_part_tmp
      REAL*8            time_d_tmp, E_part_tmp, part_dir_tmp(3)
c      REAL*8           dir, dir_tmp

      INTEGER           ith_event, jth_particle, kth_dir

      INTEGER           read_count       ! only need to read titles once

      SAVE              first_call, read_count

      DATA              first_call /.true./,
     +                  read_count /0/

**********************************************************************

**********************************************************************
*  Initializations
**********************************************************************

      DO 1 kth_dir = 1, 3
         part_dir(kth_dir) = 0.d0
         part_dir_tmp(kth_dir) = 0.d0
1     CONTINUE

      id_part_tmp = 0
      time_d_tmp  = 0.d0
      E_part_tmp  = 0.d0
      decay_delay = 0.d0

**********************************************************************

      WRITE (*,5) interaction, n_events
5     FORMAT (1x,'  Sampling interaction',i3,' to produce final state pr
     +oducts for',i7,' events')  

**************************************************************************
*  If not, then read sampled neutrino event information 
**************************************************************************

         CALL file_open(91,nu_sample_file,interaction)

         DO 10 ith_event = 1, n_events
            READ (91,*) time_d(ith_event), E_nu(ith_event)
10       CONTINUE

         CLOSE (unit = 91)

***************************************************************************
*  Send information down the correct interaction class pathway
***************************************************************************

         CALL file_open(93,part_sample_file,interaction)

         CALL get_int_type(interaction,int_type)
         CALL get_num_part(interaction,n_particles,id_part)

         n_particles_total = n_events * n_particles

* get supernova direction cosines

         IF (first_call) CALL sn_direction(sn_dir_cos)

****************************************************************************
*  Sample elastic scattering reactions event by event
*  (don't need to worry about 'id_part' because there's only one final
*  state particle, namely the scattered electron)
****************************************************************************

         IF (int_type .EQ. id_nu_es) THEN

            IF (    interaction .EQ. 58 
     +         .OR. interaction .EQ. 62
     +         .OR. interaction .EQ. 67
     +         .OR. interaction .EQ. 70
     +         .OR. interaction .EQ. 73
     +         .OR. interaction .EQ. 76) THEN

               DO 15 ith_event = 1, n_events
                  CALL en_nu_p_es(interaction,E_nu(ith_event),E_part)
                  CALL dir_nu_p_es(E_nu(ith_event),E_part,
     +                          sn_dir_cos,part_dir)

                  WRITE (93,18) id_part(1),time_d(ith_event),E_part,
     +                       part_dir(1),part_dir(2),part_dir(3)
18                FORMAT (i2,f15.8,f14.8,f14.8,f14.8,f14.8)
15             CONTINUE

            ELSE

              DO 20 ith_event = 1, n_events
                 CALL en_nu_e_es(interaction,E_nu(ith_event),E_part)
                 CALL dir_nu_e_es(E_nu(ith_event),E_part,
     +                          sn_dir_cos,part_dir)

                 WRITE (93,25) id_part(1),time_d(ith_event),E_part, 
     +                       part_dir(1),part_dir(2),part_dir(3)
25               FORMAT (i2,f15.8,f14.8,f14.8,f14.8,f14.8)
20            CONTINUE

            ENDIF

****************************************************************************
*  Sample charge current reactions event by event over each final state
*  detected particle
*  (also deal with nu_e_bar + d --> n + n + e+ kinematics)
****************************************************************************

         ELSEIF (int_type .EQ. id_nu_cc) THEN

*** Read in titles files for 'nu + d' double differential cross sections (once)
            IF (     interaction .EQ. 3
     +          .OR. interaction .EQ. 15) THEN

               IF (read_count .EQ. 0) THEN
                  CALL read_titles()
                  read_count = 1
               ENDIF

*** Define variables needed by the SNOMAN routine, 'f_xs_diff.for'
               CALL get_nu_type(interaction,nu_type)
               IQ(LMCTKP+KMCTK_IDP) = nu_type
               C_SCTN = 0 ! KN cross sections
c               C_SCTN = 5 ! NSAPMGK cross sections (they extend to 170 MeV)

            ENDIF

****************************************************************************

            DO 30 ith_event = 1, n_events
32             CONTINUE
               DO 31 jth_particle = 1, n_particles
                  CALL part_cc(interaction,E_nu(ith_event),
     +                         id_part(jth_particle),
     +                         E_part,sn_dir_cos,part_dir,redo,
     +                         decay_delay)

* check nu_e_bar + d --> n + n + e+ kinematics
                  IF (interaction .EQ. 15) THEN
* wait to see if the kinematics work for this event before writing to file
* (hold first event in buffer until we can check the first neutron)
                     IF (jth_particle .EQ. 1) THEN
                        id_part_tmp     = id_part(jth_particle)
                        time_d_tmp      = time_d(ith_event)
                        E_part_tmp      = E_part
                        part_dir_tmp(1) = part_dir(1)
                        part_dir_tmp(2) = part_dir(2)
                        part_dir_tmp(3) = part_dir(3)
c                        WRITE (*,*) 'E_e = ', E_part_tmp
                     ELSEIF (jth_particle .EQ. 2) THEN
* if the kinematics work out, then write the information for both the
* positron and the first neutron to file
                        IF (redo .EQ. 0) THEN
c                           dir_tmp =  part_dir_tmp(1)**2
c     +                              + part_dir_tmp(2)**2
c     +                              + part_dir_tmp(3)**2
c
c                           WRITE (*,36) dir_tmp
                           WRITE (93,35) id_part_tmp, 
     +                                   time_d_tmp,E_part_tmp,
     +                                   part_dir_tmp(1),
     +                                   part_dir_tmp(2),
     +                                   part_dir_tmp(3)

c                           dir =  part_dir(1)**2
c     +                          + part_dir(2)**2
c     +                          + part_dir(3)**2
c
c                           WRITE (*,36) dir

                           WRITE (93,35) id_part(jth_particle),
     +                                   time_d(ith_event),E_part, 
     +                                   part_dir(1),
     +                                   part_dir(2),
     +                                   part_dir(3)
* have to redo calculations, so go back to the beginning of this event
                        ELSE
                           WRITE (*,*) 'redoing nu_e_bar + d sample (E_e
     += ', E_part_tmp,', time_d = ', time_d_tmp,')' 
                           GOTO 32       ! I know, it's ugly :(
                        ENDIF
                     ELSE
c                        dir =  part_dir(1)**2  
c     +                       + part_dir(2)**2
c     +                       + part_dir(3)**2
c                        
c                        WRITE (*,36) dir

                        WRITE (93,35) id_part(jth_particle),
     +                                time_d(ith_event),E_part,
     +                                part_dir(1),
     +                                part_dir(2),
     +                                part_dir(3)
                     ENDIF 
                  ELSE
                     WRITE (93,35) id_part(jth_particle),
     +                             time_d(ith_event)+decay_delay,
     +                             E_part, 
     +                             part_dir(1),part_dir(2),part_dir(3)
                  ENDIF
35                FORMAT (i2,f15.8,f14.8,f14.8,f14.8,f14.8)
36                FORMAT ('direction magnitude = ', f8.5)
31             CONTINUE
30          CONTINUE

****************************************************************************
*  Sample neutral current reactions event by event
****************************************************************************

         ELSEIF (int_type .EQ. id_nu_nc) THEN

*** Read in titles files for 'nu + d' double differential cross sections (once)
            IF (     interaction .EQ. 8
     +          .OR. interaction .EQ. 19
     +          .OR. interaction .EQ. 23
     +          .OR. interaction .EQ. 32
     +          .OR. interaction .EQ. 41
     +          .OR. interaction .EQ. 50) THEN

               IF (read_count .EQ. 0) THEN
                  CALL read_titles()
                  read_count = 1
               ENDIF

*** Define variables needed by the SNOMAN routine, 'f_xs_diff.for'
               CALL get_nu_type(interaction,nu_type)
               IQ(LMCTKP+KMCTK_IDP) = nu_type
               C_SCTN = 0 ! KN cross sections
                          ! NSAPMGK cross sections do not extend past 50 MeV

            ENDIF
 
            DO 40 ith_event = 1, n_events
               DO 41 jth_particle = 1, n_particles
                  CALL part_nc(interaction,E_nu(ith_event),
     +                         id_part(jth_particle),
     +                         E_part,sn_dir_cos,part_dir)
                    WRITE (93,45) id_part(jth_particle),
     +                          time_d(ith_event),E_part,
     +                          part_dir(1),part_dir(2), part_dir(3)
45                  FORMAT (i2,f15.8,f14.8,f14.8,f14.8,f14.8)
41             CONTINUE
40          CONTINUE

         ENDIF

         CLOSE (unit = 93)

      RETURN
      END
