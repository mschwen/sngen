***********************************************************************
*  Filename: create_mcpl.f
***********************************************************************
*
*  Creates MCPL file which containing the sampled SN neutrino events
*  (not the secondary particles).  This file is read in by SNOMAN
*
*  Inputs:
*	o  interaction id number (not array)
*	o  count total for the specified interaction (not array)
*  Outputs:
*	o  MCPL file with sampled time (ns), energy (MeV)
*
*************************************************************************

      SUBROUTINE create_mcpl(interaction,num_counts)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'sample_count_size.inc'
      INCLUDE 'convert_units.inc'

C  Common Block Declarations

      COMMON /sn_interaction_names/ interaction_name
      CHARACTER*100 interaction_name(interaction_total)

      COMMON /part_sample_file/ part_sample_file
      CHARACTER*100 part_sample_file

      COMMON /mcpl_stuff/ mcpl_flag, mcpl_file, mcpl_default
      INTEGER mcpl_flag
      CHARACTER*100 mcpl_file, mcpl_default

C  Global Variable Declarations (in)

      INTEGER interaction, num_counts

C  Local Variable Declarations

      INTEGER num_part, id_part
      REAL*8 time, energy, u, v, w
      REAL*8 x, y, z

      INTEGER i

************************************************************************

************************************************************************
*  Initializations
************************************************************************

      num_part = 1   ! only 1 neutrino per MCPL entry

**** set position coordinates to 0 because we rely on SNOMAN to
**** randomly position the event within the specified region

      x = 0.
      y = 0.
      z = 0.

*************************************************************************
*  Output MCPL data list
*************************************************************************

      CALL file_open(98,mcpl_file,interaction)

C      WRITE (*,*) interaction, num_counts

*****  MCPL header

      WRITE (98,10) 
10    FORMAT ('*LOG')
      IF (interaction .LT. 10) THEN
         WRITE (98,20) interaction
20       FORMAT ('*---  Supernova Titles File: mcpl0',i1,'.out
     + 10 December 2002')
      ELSEIF (interaction .GE. 10) THEN
         WRITE (98,21) interaction
21       FORMAT ('*---  Supernova Titles File: mcpl',i2,'.out 
     + 10 December 2002')
      ENDIF

      WRITE (98,30) 
30    FORMAT ('*---                 Banks: MCPL')
      WRITE (98,40) 
40    FORMAT ('*LOG OFF')
      WRITE (98,50)
* change the last number to allow for more MCPL entries 
* (currently large enough for 100 supernovae at 10 kpc)
50    FORMAT ('*DO  MCPL   1  -i(30I 2I / 2i 8f ) -N600000')
      WRITE (98,*)
      WRITE (98,51)
51    FORMAT ('#.  This bank contains specific') 
      WRITE (98,60) 
60    FORMAT ('#.  seed vertex information which is reproduced')
      WRITE (98,70) 
70    FORMAT ('#.  verbatim by the Monte Carlo.')
      WRITE (98,*)
      WRITE (98,75) 
75    FORMAT ('#.')
      WRITE (98,80) 
80    FORMAT ('#.  MCPL file for the following supernova interaction:') 
      WRITE (98,90) interaction_name(interaction)
90    FORMAT ('#.  ', a)
      WRITE (98,100) 
100   FORMAT ('#.')
      WRITE (98,110) 
110   FORMAT ('#.    Contact: Jaret Heise,  LANL')
      WRITE (98,120) 
120   FORMAT ('#.')
      WRITE (98,130) 
130   FORMAT ('#.       Standard Database Header')
      WRITE (98,140) 
140   FORMAT ('#.')
      WRITE (98,150) 
150   FORMAT ('19750101        0 20380517 03331900  #.  1..4 Intrinsic v
     +alidity')
      WRITE (98,160) 
160   FORMAT ('       0        0        0           #.  5..7 Data type,
     +Task type, Format no.')
      WRITE (98,170) 
170   FORMAT ('       0        0        0           #.  8..10 Creation D
     +ate, Time, Source Id.')
      WRITE (98,180) 
180   FORMAT ('19750101        0 20380517 03331900  #. 11..14 Effective
     +validity')
      WRITE (98,190) 
190   FORMAT ('       0        0                    #. 15..16 Entry Date
     +Time')
      WRITE (98,200) 
200   FORMAT ('4*0                                  #. 17..20 Spare')
      WRITE (98,210) 
210   FORMAT ('10*0                                 #. 21..30 Temporary
     +data (not in database)')
      WRITE (98,220) 
220   FORMAT ('#.')
      WRITE (98,230) 
230   FORMAT ('#.    End of Standard Database Header')
      WRITE (98,240) 
240   FORMAT ('#.')
      WRITE (98,*)
      WRITE (98,250) 
250   FORMAT ('#.    User Data.')

***** MCPL Data

* number of events with 10 words per event
      WRITE (98,260) num_counts
260   FORMAT (i6,' 10')
      WRITE (98,*)
      WRITE (98,265)
265   FORMAT ('#. id  x   y   z  time (ns)         E (MeV)     dircos_x
     +   dircos_y    dircos_z')

***** read in sampled final state particle data and output MCPL data
***** (remember to convert time to "ns" from "ms" for SNOMAN)

      CALL file_open(97,part_sample_file,interaction)

      DO 270 i = 1, num_counts
         READ (97,*) id_part, time, energy, u, v, w
         IF (energy .NE. 0) THEN
          WRITE (98,275) num_part, id_part, x, y, z, 
     +                  (time * ms2ns), energy, u, v, w
275       FORMAT (i1,2x,i2,3(1x,f3.0),1x,e17.10,1x,f12.8,3(1x,f11.8))
         ENDIF
270   CONTINUE

      CLOSE (unit = 97)

*************************************************************************

      CLOSE (unit = 98)

      RETURN
      END
