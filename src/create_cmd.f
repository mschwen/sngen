***********************************************************************
*  Filename: create_cmd.f
***********************************************************************
*
*  Creates .cmd files which are used by SNOMAN.  Current settings include:
*       o no cross talk
*       o measured trigger efficiency curve, DQXX banks
*
*  Inputs:
*	o  supernova direction cosines
*	o  interaction id number (not array)
*	o  count total for the specified interaction (not array)
*
*************************************************************************

      SUBROUTINE create_cmd(interaction,num_counts)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'id_part.inc'
      INCLUDE 'id_cmd_type.inc'

C  Common Block Declarations

      COMMON           /sn_interaction_names/
     +                 interaction_name
      CHARACTER*57     interaction_name(interaction_total)

      COMMON           /mcpl_stuff/
     +                 mcpl_flag, mcpl_file, mcpl_default
      INTEGER          mcpl_flag
      CHARACTER*40     mcpl_file, mcpl_default

      COMMON           /cmd_stuff/ 
     +                 cmd_flag, cmd_default,
     +                 zdab_flag, zdab_default,
     +                 hbook_default
      INTEGER          cmd_flag, zdab_flag
      CHARACTER*40     cmd_default, zdab_default, hbook_default

      COMMON           /snoman_parameters/
     +                 nhit_thres, comp_thres,
     +                 pmt_efficiency,
     +                 h2o_in_d2o, O17_in_d2o
      REAL*8           nhit_thres, comp_thres,
     +                 pmt_efficiency,
     +                 h2o_in_d2o, O17_in_d2o

C  Global Variable Declarations

      INTEGER          interaction, num_counts

C  Local Variable Declarations

      CHARACTER*40     medium

      CHARACTER*40     cmd_int_file
      CHARACTER*40     mcpl_int_file
      CHARACTER*40     zdab_int_file
      CHARACTER*40     hbook_int_file

      INTEGER          i, j
      INTEGER          salt_flag, cmd_types, type
      INTEGER          n_part, id_part(3)

************************************************************************
*  Initializations
************************************************************************

      CALL get_medium(interaction,medium)

************************************************************************
*  Check to see how many cmd files are required for this interaction
*  (ie., whether we need to specify "nosalt", "salt" and "ncd" files;
*  since the NCD strings interfere with Cerenkov light propagation,
*  all interactions will require "ncd" files; not all interactions
*  require "salt")
************************************************************************

      salt_flag = 0  ! only some interactions will require a salt cmd file

      CALL get_num_part(interaction,n_part,id_part)
     
      DO 4 i = 1, n_part
        IF (id_part(i) .EQ. id_neutron) THEN
           salt_flag = 1
        ENDIF
4     CONTINUE

      IF (salt_flag .GE. 1) THEN
        cmd_types = 3  ! pure D2O, salt, NCDs
      ELSE
        cmd_types = 2  ! pure D2O, NCDs
      ENDIF

**************************************************************************
*  Loop over the number of types (nosalt/salt, ncd or no distinction) 
*  of cmd output files required
**************************************************************************

      DO 5 j = 1, cmd_types

         CALL file_increment(mcpl_default,interaction,mcpl_int_file)
         CALL file_increment(zdab_default,interaction,zdab_int_file)
         CALL file_increment(hbook_default,interaction,hbook_int_file)

         type = (j - 1)  ! 0=pure D2O, 1=NCDs, 2=salt

         IF (      type .EQ. 0              ! pure D2O
     +       .AND. cmd_types .LT. 3) THEN   ! no neutrons
            CALL file_open(99,cmd_default,interaction) ! just use generic
         ELSE
            CALL file_increment(cmd_default,interaction,cmd_int_file)

            CALL file_cmd_type(cmd_int_file,type,cmd_int_file)
            CALL file_cmd_type(zdab_int_file,type,zdab_int_file)
            CALL file_cmd_type(hbook_int_file,type,hbook_int_file)

            OPEN (unit = 99, file = cmd_int_file, status = 'unknown')
         ENDIF

*************************************************************************
*  Output .cmd commands
*************************************************************************

      WRITE (99,10) 
10    FORMAT ('*********************************************************
     +************')

      IF (interaction .LT. 10) THEN
         WRITE (99,11) interaction
11       FORMAT ('*  Supernova Command File: sn0',i1,'.cmd
     +   18 March 2003')
      ELSEIF (interaction .GE. 10) THEN
         WRITE (99,12) interaction
12       FORMAT ('*  Supernova Command File: sn',i2,'.cmd
     +   18 March 2003')
      ENDIF

      WRITE (99,13)
13    FORMAT ('*')
      WRITE (99,14)
14    FORMAT ('*  Contact: Jaret Heise, LANL')
      WRITE (99,15)
15    FORMAT ('*')
      WRITE (99,20) 
20    FORMAT ('*  SNOMAN cmd file for the supernova interaction:')
      WRITE (99,30) interaction_name(interaction)
30    FORMAT ('*  ', a)
      WRITE (99,40) 
40    FORMAT ('*********************************************************
     +************')

      WRITE (99,*)

      WRITE (99,50)
50    FORMAT ('***  Job Control')
      WRITE (99,*)
      WRITE (99,60)
60    FORMAT ('$processor_list ''MCO UCL PCK OUT END''')
      WRITE (99,*)
      WRITE (99,80)
80    FORMAT ('$old_eca')
      WRITE (99,100)
100   FORMAT ('set bank tcal 1 word 4 to 1')
      WRITE (99,109)
109   FORMAT ('titles cal_const1_mc.dat')
      WRITE (99,110)
110   FORMAT ('titles ccc_info')
      WRITE (99,*)
      IF (zdab_flag .EQ. 1) THEN
         WRITE (99,141) zdab_int_file
141      FORMAT ('file out 1 ', a,'format=$zdab_only')
         WRITE (99,142)
142      FORMAT ('$zdab_option $zdab_min_mc')
         WRITE (99,*)
         WRITE (99,143)
143      FORMAT ('*$disable_ntuple 500')
         WRITE (99,144)
144      FORMAT ('*$enable_ntuple  600')
         WRITE (99,145) 
145      FORMAT ('*titles my_ntuple.dat')
         WRITE (99,*)
         WRITE (99,146) hbook_int_file
146      FORMAT ('*$hbook_file            ''', a,'''')
         WRITE (99,*)
         WRITE (99,147)
147      FORMAT ('$seed_file             ''rseed.dat        ''')
         WRITE (99,148)
148      FORMAT ('$pegs_file             ''pegs4_10.dat     ''')
      ELSE
         WRITE (99,151) zdab_int_file
151      FORMAT ('*file out 1 ', a,'format=$zdab_only')
         WRITE (99,152)                          
152      FORMAT ('*$zdab_option $zdab_min_mc')
         WRITE (99,*)
         WRITE (99,153)
153      FORMAT ('$disable_ntuple 500')
         WRITE (99,154)
154      FORMAT ('$enable_ntuple  600')
         WRITE (99,155) 
155      FORMAT ('titles my_ntuple.dat')
         WRITE (99,*)
         WRITE (99,156) hbook_int_file
156      FORMAT ('$hbook_file            ''', a,'''')
         WRITE (99,*)
         WRITE (99,157)
157      FORMAT ('$seed_file             ''rseed.dat        ''')
         WRITE (99,158)
158      FORMAT ('$pegs_file             ''pegs4_10.dat     ''')
      ENDIF
      WRITE (99,*)
      WRITE (99,160)
160   FORMAT ('***  MC Control')
      WRITE (99,*)
      WRITE (99,170) interaction
170   FORMAT ('$mcrun                 ', i2,'   * label run with interac
     +tion number')
      WRITE (99,180)
180   FORMAT ('$mc_num_seed_vx        1') 
      WRITE (99,190)
190   FORMAT ('$killvx                7')
      WRITE (99,200) num_counts
200   FORMAT ('$event_limit          ', i6)
      WRITE (99,207)
207   FORMAT ('$initial_date          20000120 07160733   * Quiet (run #
     +10879_001)')
c      WRITE (99,208)
c208   FORMAT ('*$initial_date          19990512 09005380   * Shark Fins 2
c     +(run #3358_000)')
c      WRITE (99,209)
c209   FORMAT ('*$initial_date          19990511 00550020   * Flashers/Muo
c     +ns (run #3323_000)') 
c      WRITE (99,210)
c210   FORMAT ('*$initial_date          19990506 01501020   * Shark Fins
c     +(run #3178_000)')
c      WRITE (99,211)
c211   FORMAT ('*$initial_date          19981013 23211224   * Breakdown1 (
c     +run #1698_000)')
c      WRITE (99,212)
c212   FORMAT ('*$initial_date          19981014 03152437   * Breakdown2 (
c     +run #1698_001)')
c      WRITE (99,213)
c213   FORMAT ('*$initial_date          19960101 07174637   * SNO Games (
c     +fake_data_01)')

      WRITE (99,*)

      WRITE (99,230)
230   FORMAT ('***  MC Generation')
      WRITE (99,*)
      WRITE (99,231)
231   FORMAT ('* set mctype to -ve in order to use MCPL bank')
      WRITE (99,232) 
232   FORMAT ('set bank mcma 1 word 11 to -1')
      WRITE (99,*)
      WRITE (99,233) mcpl_int_file
233   FORMAT ('titles ', a)
      WRITE (99,*)
      WRITE (99,241)
241   FORMAT ('$mc_interaction_type   $pick_particle_list')
      WRITE (99,250) medium
250   FORMAT ('$mc_position           $pos_region ',a)
      WRITE (99,260) 
260   FORMAT ('$mc_time               $tim_mcpl')
      WRITE (99,270)
270   FORMAT ('$mc_energy             $en_mcpl')
      WRITE (99,290)
290   FORMAT ('$mc_direction          $dir_mcpl')

      WRITE (99,*)

      WRITE (99,480)
480   FORMAT ('***  Cerenkov Light Generation')
      WRITE (99,*)
      WRITE (99,490)
490   FORMAT ('$low_wave_cut          200.e-7')
      WRITE (99,500)
500   FORMAT ('$high_wave_cut         720.e-7')
      WRITE (99,510)
510   FORMAT ('$mode_cerfac           1      ')
      WRITE (99,511)
511   FORMAT ('$cerfac                0.28   ')
      WRITE (99,520)
520   FORMAT ('$cerfac_ii             0.35   ')

      WRITE (99,*)

      WRITE (99,530)
530   FORMAT ('***  Particle Transport')
      WRITE (99,*)
      WRITE (99,540)
540   FORMAT ('$fresnel_scat          $on')
      WRITE (99,550)
550   FORMAT ('$rayleigh_scat         $on')
      WRITE (99,560)
560   FORMAT ('$photodisintegration   $off')

      WRITE (99,*)

      WRITE (99,569)
569   FORMAT ('***  SNO Regions')
      WRITE (99,*)
      WRITE (99,570)
570   FORMAT ('$num_d2o_regions 2')
      WRITE (99,571)
571   FORMAT ('$d2o_region_01 to $D2O')
      WRITE (99,572)
572   FORMAT ('$d2o_region_02 to $ACRC_IVL')   
      WRITE (99,573) h2o_in_d2o
573   FORMAT ('$d2o_region_01 $h2o_fraction ', f8.6)
      WRITE (99,574) h2o_in_d2o
574   FORMAT ('$d2o_region_02 $h2o_fraction ', f8.6)
      IF (type .EQ. id_salt) THEN
         WRITE (99,575)
575      FORMAT ('$d2o_region_01 $salt_gmgm $added')
         WRITE (99,576)
576      FORMAT ('$d2o_region_02 $salt_gmgm $added')
      ELSE
         WRITE (99,577)
577      FORMAT ('$d2o_region_01 $salt_gmgm $nosalt')
         WRITE (99,578)
578      FORMAT ('$d2o_region_02 $salt_gmgm $nosalt')
      ENDIF
      WRITE (99,579) O17_in_d2o
579   FORMAT ('$set_O17_d2o ', f8.6)
      IF (type .EQ. id_ncd) THEN
         WRITE (99,580)
580      FORMAT ('@geom_ncd')
      ENDIF

      WRITE (99,*)

      WRITE (99,600)
600   FORMAT ('***  MC PMT Simulation / Optics')
      WRITE (99,*)
c      WRITE (99,601)
c601   FORMAT ('$grey_disk_pmt_bounce  $on')
c      WRITE (99,602)
c602   FORMAT ('$grey_disk_pmt')
      WRITE (99,603)
603   FORMAT ('$3d_pmt')
      WRITE (99,604)
604   FORMAT ('* Enable modified 3D-PMT angular response')
      WRITE (99,605)
605   FORMAT ('set bank GEDS 423 word 7  to   1 * Select optics6b')
      WRITE (99,606)
606   FORMAT ('set bank GEDS 423 word 28 to 0.4 * set stack reflectivity
     +')
      WRITE (99,607)
607   FORMAT ('set bank GEDS 423 word 29 to 0.725 *set pcath thick reduc
     +. factor')
      WRITE (99,608)
608   FORMAT ('$acrylic_var $off')
      WRITE (99,609)
609   FORMAT ('@geom_belly_plate_rope')
      WRITE (99,610)
610   FORMAT ('@geom_acrylic_tiles')
      WRITE (99,611)
611   FORMAT ('* Use Bryce''s latest media.dat settings')
      WRITE (99,612)
612   FORMAT ('@media_qoca_jan22')
c      WRITE (99,613)
c613   FORMAT ('* Load greydisk PMT angular response')
c      WRITE (99,614)
c614   FORMAT ('titles pmt_response_qoca_jan22.dat')
      WRITE (99,615) pmt_efficiency
615   FORMAT ('set bank MCMA 1 word 24 to ', f6.4)
      WRITE (99,616)
616   FORMAT ('* PMT noise rate for Run #10879')
      WRITE (99,617)
617   FORMAT ('$pmt_nois_rate 568.01')

      WRITE (99,*)

      WRITE (99,620)
620   FORMAT ('***  MC DAQ')
      WRITE (99,*)
      WRITE (99,621) nhit_thres
621   FORMAT ('$mcdaq_hthresh        ', f5.1)
      WRITE (99,622) comp_thres
622   FORMAT ('$mcdaq_cthresh        ', f5.1)
      WRITE (99,623)
623   FORMAT ('$mcdaq_ztell           $on')
      WRITE (99,624)
624   FORMAT ('*$mcdaq_xtalk_thresh    30.')
      WRITE (99,625)
625   FORMAT ('$trigger_mode          $simulated_dead')

      WRITE (99,*)

      WRITE (99,670)
670   FORMAT ('***  Job Specification Complete: read titles files and ex
     +ecute')
      WRITE (99,*)
c      WRITE (99,675)
c675   FORMAT ('@read_neutron_titles')
      WRITE (99,680) 
680   FORMAT ('@run_snodb')
      WRITE (99,690)
690   FORMAT ('*$go')

      CLOSE (unit = 99)

5     CONTINUE

*************************************************************************

      RETURN
      END
