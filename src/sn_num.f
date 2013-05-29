***********************************************************************
*  Filename: sn_num.f
***********************************************************************
*
*  Inputs:
*	o  nu_type
*	o  neutrino time (ms)
*  Outputs:
*       o  number of neutrinos [L_nu(t) / <E_nu>]
*
*  Note: this routine only depends on time, so when integrating, we should
*  save the previous number 
* 
***********************************************************************
     
      FUNCTION sn_num(nu_type,this_time)

      IMPLICIT none

      INCLUDE 'id_nu.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON       /sn_lumen_fit_file/
     +             lumen_fit_file_flag,
     +             lumen_nu_e_file, lumen_nu_eb_file, lumen_nu_mu_file
      INTEGER      lumen_fit_file_flag
      CHARACTER*40 lumen_nu_e_file, lumen_nu_eb_file, lumen_nu_mu_file

      COMMON       /sn_time_entries/
     +             time_entries1, time_entries2, time_entries
      INTEGER      time_entries1, time_entries2, time_entries 

C  Global Variable Declarations

      INTEGER      interaction            ! input
      REAL*8       this_time              ! input

C  Local Variable Declarations

      REAL*8       sn_num, sn_lum, sn_Eavg,    ! functions
     +             lum, Eavg

      LOGICAL      done_nu_e, done_nu_eb, done_nu_mu

      REAL*8       last_time, last_num

      INTEGER      n_calls, last_interaction, last_nu, nu_type

      SAVE         n_calls,
     +             last_interaction, last_nu, last_time, last_num,
     +             done_nu_e, done_nu_eb, done_nu_mu

************************************************************************
*  Initializations
************************************************************************

      DATA         n_calls          /0/,
     +             last_interaction /0/,
     +             last_time        /-99999./,
     +             last_nu          /0/,
     +             last_num         /0./

      DATA         done_nu_e        /.false./,  ! luminosity/<E> output
     +             done_nu_eb       /.false./,  ! luminosity/<E> output
     +             done_nu_mu       /.false./   ! luminosity/<E> output

**************************************************************************

C      CALL get_nu_type(interaction,nu_type)

      IF (      this_time .EQ. last_time
     +    .AND. nu_type   .EQ. last_nu  ) THEN

         sn_num = last_num  ! if same time and neutrino type, use last sn_num

      ELSE   ! calculate sn_num

*************************************************************************
*  Calculate number of supernova neutrinos at 'time' assuming an energy
*  distribution centred around E_nu_ave
*    o sn_lum  --> supernova neutrino luminosity (ergs/s)
*    o sn_Eavg --> supernova neutrino average energy (MeV)
*************************************************************************

         lum = sn_lum(nu_type,this_time)
         Eavg = sn_Eavg(nu_type,this_time)

         sn_num = (erg2mev * lum) / Eavg

**************************************************************************
*  Write the various luminosity and energy spectra to file.  
*  -> gets here if this_time .NE. last_time; ie., skip energy integration
*  (assumes time integration in "count_sum" done first, then energy)
**************************************************************************
      
         IF (lumen_fit_file_flag .EQ. 1) THEN

* reset counter for the next interaction

            IF (interaction .NE. last_interaction) n_calls = 0

            n_calls = n_calls + 1       ! increment time step counter

            IF (.NOT. done_nu_e
     +          .AND. nu_type .EQ. id_nu_e) THEN

               WRITE (60,*) this_time, lum, Eavg

* don't need to write these values again after we've written the
* entire set out once (same goes for others)
               
               IF (n_calls .EQ. time_entries) done_nu_e = .true.

            ELSEIF (.NOT. done_nu_eb
     +              .AND. nu_type .EQ. id_nu_ebar) THEN

               WRITE (61,*) this_time, lum, Eavg
               
               IF (n_calls .EQ. time_entries) done_nu_eb = .true.

            ELSEIF (.NOT. done_nu_mu
     +              .AND. nu_type .EQ. id_nu_mu) THEN

               WRITE (62,*) this_time, lum, Eavg
               
               IF (n_calls .EQ. time_entries) done_nu_mu = .true.
            ENDIF ! check for neutrino type
         ENDIF    ! check for file_flag

* remember a few things about the last time this routine was called

         last_interaction = interaction
         last_nu = nu_type
         last_time = this_time
         last_num = sn_num

      ENDIF  ! use last sn_num value

      END
