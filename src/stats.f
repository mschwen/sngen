***************************************************************************
*  Filename: stats.f
***************************************************************************

      SUBROUTINE stats(exp_id,id_interaction,count_pois,icount_pois)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'id_nu.inc'
      INCLUDE 'time_size.inc'
      INCLUDE 'sample_count_size.inc'

C  Global Variable Declarations

      INTEGER          id_interaction(*),     ! input
     +                 icount_pois(*),        ! input
     +                 exp_id                 ! input
      REAL*8           count_pois(*)          ! input

C  Common Block Declarations

      COMMON           /sn_model_stats/
     +                 model_stats
      INTEGER          model_stats

      COMMON           /generator_flag/
     +                 gen_flag, counts_file
      INTEGER          gen_flag
      CHARACTER*40     counts_file

      COMMON           /sample_stuff/
     +                 sample_E_nu
      REAL*8           sample_E_nu(interaction_total,sample_count_size)

      COMMON           /stats_stuff/
     +                 t_d, cr_d, c_d, cumc_d

      REAL*8           t_d(time_d_size),
     +                 cr_d(interaction_total,time_d_size),
     +                 c_d(interaction_total,time_d_size),
     +                 cumc_d(interaction_total,time_d_size)
 
      COMMON           /nu_countrate_file/
     +                 all_crfile,
     +                 es_crfile, cc_crfile, nc_crfile,
     +                 nu_e_crfile, nu_eb_crfile, nu_mu_crfile,
     +                 d2o_crfile, h2o_crfile, av_crfile
      CHARACTER*40     all_crfile,
     +                 es_crfile, cc_crfile, nc_crfile,
     +                 nu_e_crfile, nu_eb_crfile, nu_mu_crfile,
     +                 d2o_crfile, h2o_crfile, av_crfile

      COMMON           /nu_count_file/
     +                 all_cfile,
     +                 es_cfile, cc_cfile, nc_cfile,
     +                 nu_e_cfile, nu_eb_cfile, nu_mu_cfile,
     +                 d2o_cfile, h2o_cfile, av_cfile
      CHARACTER*40     all_cfile,
     +                 es_cfile, cc_cfile, nc_cfile,
     +                 nu_e_cfile, nu_eb_cfile, nu_mu_cfile,
     +                 d2o_cfile, h2o_cfile, av_cfile

      COMMON           /nu_cumulativecount_file/
     +                 all_cumcfile,
     +                 es_cumcfile, cc_cumcfile, nc_cumcfile,
     +                 nu_e_cumcfile, nu_eb_cumcfile, nu_mu_cumcfile,
     +                 d2o_cumcfile, h2o_cumcfile, av_cumcfile
      CHARACTER*40     all_cumcfile,
     +                 es_cumcfile, cc_cumcfile, nc_cumcfile,
     +                 nu_e_cumcfile, nu_eb_cumcfile, nu_mu_cumcfile,
     +                 d2o_cumcfile, h2o_cumcfile, av_cumcfile

      COMMON           /interaction_group/
     +                 pick_group_es, pick_group_cc, pick_group_nc,
     +                 pick_group_nu_e,
     +                 pick_group_nu_eb,
     +                 pick_group_nu_mu,
     +                 pick_group_d2o, pick_group_h2o, pick_group_av
      INTEGER          pick_group_es, pick_group_cc, pick_group_nc,
     +                 pick_group_nu_e,
     +                 pick_group_nu_eb,
     +                 pick_group_nu_mu,
     +                 pick_group_d2o, pick_group_h2o, pick_group_av

      COMMON           /interaction_arrays/
     +                 interaction_all,
     +                 interaction_es, interaction_cc, interaction_nc,
     +                 interaction_nu_e, interaction_nu_eb, 
     +                 interaction_nu_mu,
     +                 interaction_d2o, interaction_h2o, 
     +                 interaction_av, interaction_LAB

      INTEGER          interaction_all(interaction_total)

      INTEGER          interaction_es(interaction_total),
     +                 interaction_cc(interaction_total),
     +                 interaction_nc(interaction_total)

      INTEGER          interaction_nu_e(interaction_total),
     +                 interaction_nu_eb(interaction_total),
     +                 interaction_nu_mu(interaction_total)

      INTEGER          interaction_d2o(interaction_total),
     +                 interaction_h2o(interaction_total),
     +                 interaction_av(interaction_total),
     +                 interaction_lab(interaction_total)

C  Local Variable Declarations

      REAL*8           e_E, e_num, e_Eavg,
     +                 ebar_E, ebar_num, ebar_Eavg,
     +                 muon_E, muon_num, muon_Eavg,
     +                 muonbar_E, muonbar_num, muonbar_Eavg
      REAL*8           muon_total_E, muon_total_num, muon_total_Eavg
      REAL*8           sn_total_E, sn_total_num, sn_total_Eavg

      REAL*8           sum_E_nu
      INTEGER          i,
     +                 inum_nu,
     +                 id, this_id(interaction_total)

      INTEGER          izero, iunit, iten
      CHARACTER*2      cid

      CHARACTER*40     this_crfile, this_cfile, this_cumcfile

      DATA             muon_total_E    /0.d0/,
     +                 muon_total_num  /0.d0/,
     +                 muon_total_Eavg /0.d0/
      DATA             sn_total_E      /0.d0/, 
     +                 sn_total_num    /0.d0/,
     +                 sn_total_Eavg   /0.d0/

**************************************************************************

* check if stats flag was set

      IF (model_stats .EQ. 0) RETURN

**************************************************************************
*  Initializations
**************************************************************************

      DO i = 1, interaction_total
         this_id(i) = 0
      END DO

      OPEN (unit = 15, file = 'stats.out', status = 'unknown')

**************************************************************************
*  Total energies and total average energies
**************************************************************************

      WRITE (15,*) 'Supernova Source Stats'
      WRITE (15,*) '~~~~~~~~~~~~~~~~~~~~~~'

* nu_e
      CALL stats_src(id_nu_e,e_E,e_num,e_Eavg)
      WRITE (15,*) 'Total electron neutrino energy is ',
     +             e_E,' ergs'
      WRITE (15,*) 'Total number of electron neutrinos is ',
     +             e_num
      WRITE (15,*) 'Total average electron neutrino energy is ',
     +             e_Eavg,' MeV'

* nu_ebar
      CALL stats_src(id_nu_ebar,ebar_E,ebar_num,ebar_Eavg)
      WRITE (15,*) 'Total electron antineutrino energy is ',
     +             ebar_E,' ergs'
      WRITE (15,*) 'Total number of electron antineutrinos is ',
     +             ebar_num
      WRITE (15,*) 'Total average electron antineutrino energy is ', 
     +             ebar_Eavg,' MeV'

* nu_mu/tau
      CALL stats_src(id_nu_mu,muon_E,muon_num,muon_Eavg) 
      WRITE (15,*) 'Total muon/tau neutrino energy is ',
     +             muon_E,' ergs'
      WRITE (15,*) 'Total number of muon/tau neutrinos is ',
     +             muon_num
      WRITE (15,*) 'Total average muon/tau neutrino energy is ', 
     +             muon_Eavg,' MeV'

* nu_mubar/taubar
      CALL stats_src(id_nu_mubar,muonbar_E,muonbar_num,muonbar_Eavg) 
      WRITE (15,*) 'Total muon/tau antineutrino energy is ',
     +             muonbar_E,' ergs'
      WRITE (15,*) 'Total number of muon/tau antineutrinos is ',
     +             muonbar_num
      WRITE (15,*) 'Total average muon/tau antineutrino energy is ', 
     +             muonbar_Eavg,' MeV'

* "nu_mu"
      muon_total_E = muon_E + muonbar_E
      muon_total_num = muon_num + muonbar_num
      muon_total_Eavg = (muon_num * muon_Eavg
     +                 + muonbar_num * muonbar_Eavg)
     +                 / (muon_num + muonbar_num)

      WRITE (15,*) 'Total "muon" neutrino energy is ',
     +             muon_total_E,' ergs'
      WRITE (15,*) 'Total number of "muon" neutrinos is ',
     +             muon_total_num
      WRITE (15,*) 'Total average "muon" neutrino energy is ',
     +             muon_total_Eavg,' MeV'

      WRITE (15,*)

      sn_total_E = e_E + ebar_E + muon_E + muonbar_E
      sn_total_num = e_num + ebar_num + muon_num + muonbar_num
      sn_total_Eavg = (e_num * e_Eavg
     +               + ebar_num * ebar_Eavg
     +               + muon_num * muon_Eavg
     +               + muonbar_num * muonbar_Eavg)
     +               / sn_total_num

      WRITE (15,*) 'Total supernova energy is ', sn_total_E,' ergs'
      WRITE (15,*) 'Total number of supernova neutrinos is ',
     +              sn_total_num
      WRITE (15,*) 'Total supernova average neutrino energy is ',
     +              sn_total_Eavg,' MeV'

*************************************************************************
*  Sample Count Stats
*************************************************************************

      IF (gen_flag .GT. 1) THEN  ! we haven't properly sampled
         WRITE (*,*)
         WRITE (*,*) 'Skipping sample stats...'

      ELSE
**************************************************************************
*  Sampled neutrino energies
**************************************************************************

      WRITE (15,*)
      WRITE (15,*) 'Detected Supernova Neutrino Stats'
      WRITE (15,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

* all interactions

      CALL group_sum2d(interaction_all,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_all,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of all sampled neutrinos is ',
     +             sum_E_nu / float(inum_nu),' MeV'

* nu_e interactions

      CALL group_sum2d(interaction_nu_e,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_nu_e,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled electron neutrinos i
     +s ', sum_E_nu / float(inum_nu),' MeV'

* nu_ebar interactions

      CALL group_sum2d(interaction_nu_eb,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_nu_eb,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled electron antineutrin
     +os is ', sum_E_nu / float(inum_nu),' MeV'

* nu_"mu" interactions

      CALL group_sum2d(interaction_nu_mu,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_nu_mu,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled "muon" neutrinos is
     + ', sum_E_nu / float(inum_nu),' MeV'

* CC interactions
      CALL group_sum2d(interaction_cc,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_cc,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled CC neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

* NC interactions

      CALL group_sum2d(interaction_nc,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_nc,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled NC neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

* ES interactions

      CALL group_sum2d(interaction_es,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_es,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled ES neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

      IF (exp_id .EQ. 1) THEN

* D2O interactionsa

         CALL group_sum2d(interaction_d2o,sample_E_nu,sum_E_nu)
         CALL group_isum1d(interaction_d2o,icount_pois,inum_nu)

         WRITE (15,*) 'Total average energy of sampled D2O neutrinos
     +   is ',sum_E_nu / float(inum_nu),' MeV'

      ELSEIF (exp_id .EQ. 2) THEN

* LAB interactions

         CALL group_sum2d(interaction_LAB,sample_E_nu,sum_E_nu)
         CALL group_isum1d(interaction_LAB,icount_pois,inum_nu)

         WRITE (15,*) 'Total average energy of sampled LAB neutrinos
     +    is ',sum_E_nu / float(inum_nu),' MeV'
 
      ENDIF

      IF (exp_id .EQ. 1 .OR. exp_id .EQ. 2) THEN

* H2O interactions

      CALL group_sum2d(interaction_h2o,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_h2o,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled H2O neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

* AV interactions

      CALL group_sum2d(interaction_av,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_av,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled AV neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

      ENDIF

      IF (exp_id .EQ. 3 ) THEN

      CALL group_sum2d(interaction_all,sample_E_nu,sum_E_nu)
      CALL group_isum1d(interaction_all,icount_pois,inum_nu)

      WRITE (15,*) 'Total average energy of sampled Pb neutrinos is ',
     + sum_E_nu / float(inum_nu),' MeV'

      ENDIF

      ENDIF  ! check for sampling

*************************************************************************
*  Count Stats
*************************************************************************

      IF (gen_flag .EQ. 3) THEN  ! we haven't filled the count arrays
         WRITE (*,*)
         WRITE (*,*) 'Skipping count stats...'
         RETURN
      ENDIF

*************************************************************************
*  Find the total count rate and cumulative counts as functions of time
*  for the specified interactions
*************************************************************************

* calculate the count stats for each individual interaction

      DO i = 1, interaction_total

         id = id_interaction(i)
         IF (count_pois(id) .GT. 0.d0) THEN   ! only do chosen reactions

* assign array element (array does not require spacing of elements)

            this_id(1) = id

* interaction ID number as a character string

            izero = ichar('0')
            iunit = mod(id,10)
            iten = mod(id-iunit,100)/10
            cid = char(izero+iten)//char(izero+iunit)

            this_crfile  = 'countrate'//cid//'.out'
            this_cfile   = 'count'//cid//'.out'
            this_cumcfile = 'cumcount'//cid//'.out'

            CALL group_sum2d_file(this_id,cr_d,t_d,this_crfile)
            CALL group_sum2d_file(this_id,c_d,t_d,this_cfile)
            CALL group_sum2d_file(this_id,cumc_d,t_d,this_cumcfile)
         ENDIF
      END DO

* find the count stats as a function of time for the group:
* all specified interactions

      CALL group_sum2d_file(interaction_all,cr_d,t_d,all_crfile)
      CALL group_sum2d_file(interaction_all,c_d,t_d,all_cfile)
      CALL group_sum2d_file(interaction_all,cumc_d,t_d,all_cumcfile)

* find the count rates and counts as a function of time for the groups:
* es, cc, nc

      IF (pick_group_es .EQ. 1) THEN
         CALL group_sum2d_file(interaction_es,cr_d,t_d,es_crfile)
         CALL group_sum2d_file(interaction_es,c_d,t_d,es_cfile)
         CALL group_sum2d_file(interaction_es,cumc_d,t_d,es_cumcfile)
      ENDIF

      IF (pick_group_cc .EQ. 1) THEN
         CALL group_sum2d_file(interaction_cc,cr_d,t_d,cc_crfile)
         CALL group_sum2d_file(interaction_cc,c_d,t_d,cc_cfile)
         CALL group_sum2d_file(interaction_cc,cumc_d,t_d,cc_cumcfile)
      ENDIF

      IF (pick_group_nc .EQ. 1) THEN
          CALL group_sum2d_file(interaction_nc,cr_d,t_d,nc_crfile)
          CALL group_sum2d_file(interaction_nc,c_d,t_d,nc_cfile)
          CALL group_sum2d_file(interaction_nc,cumc_d,t_d,nc_cumcfile)
      ENDIF

* find the count rates and counts as a function of time for the groups:
* nu_e, nu_eb, `nu_mu'

      IF (pick_group_nu_e .EQ. 1) THEN
          CALL group_sum2d_file(interaction_nu_e,cr_d,t_d,
     +                          nu_e_crfile)
          CALL group_sum2d_file(interaction_nu_e,c_d,t_d,
     +                          nu_e_cfile)
          CALL group_sum2d_file(interaction_nu_e,cumc_d,t_d,
     +                          nu_e_cumcfile)
      ENDIF

      IF (pick_group_nu_eb .EQ. 1) THEN
          CALL group_sum2d_file(interaction_nu_eb,cr_d,t_d,
     +                          nu_eb_crfile)
          CALL group_sum2d_file(interaction_nu_eb,c_d,t_d,
     +                          nu_eb_cfile)
          CALL group_sum2d_file(interaction_nu_eb,cumc_d,t_d,
     +                     nu_eb_cumcfile)
      ENDIF

      IF (pick_group_nu_mu .EQ. 1) THEN
          CALL group_sum2d_file(interaction_nu_mu,cr_d,t_d,
     +                          nu_mu_crfile)
          CALL group_sum2d_file(interaction_nu_mu,c_d,t_d,
     +                          nu_mu_cfile)
          CALL group_sum2d_file(interaction_nu_mu,cumc_d,t_d,
     +                     nu_mu_cumcfile)
      ENDIF

* find the count rates and counts as a function of time for the groups:
* D2O, H2O, AV

      IF (pick_group_d2o .EQ. 1) THEN
          CALL group_sum2d_file(interaction_d2o,cr_d,t_d,d2o_crfile)
          CALL group_sum2d_file(interaction_d2o,c_d,t_d,d2o_cfile)
          CALL group_sum2d_file(interaction_d2o,cumc_d,t_d,
     +                          d2o_cumcfile)
      ENDIF

      IF (pick_group_h2o .EQ. 1) THEN
          CALL group_sum2d_file(interaction_h2o,cr_d,t_d,h2o_crfile)
          CALL group_sum2d_file(interaction_h2o,c_d,t_d,h2o_cfile)
          CALL group_sum2d_file(interaction_h2o,cumc_d,t_d,
     +                          h2o_cumcfile)
      ENDIF

      IF (pick_group_av .EQ. 1) THEN
          CALL group_sum2d_file(interaction_av,cr_d,t_d,av_crfile)
          CALL group_sum2d_file(interaction_av,c_d,t_d,av_cfile)
          CALL group_sum2d_file(interaction_av,cumc_d,t_d,
     +                          av_cumcfile)
      ENDIF

      CLOSE (unit = 15)

      RETURN
      END
