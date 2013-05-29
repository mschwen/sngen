*************************************************************************
*  Filename: summary.f
*************************************************************************
*
*  Inputs: 
*       o array of average counts for all specified interactions
*  Outputs:
*       o summary output file
*
*************************************************************************

      SUBROUTINE summary(id_interaction,count_pois,icount_pois,
     +                   count_part,icount_part,exp_id)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'sn_classes.inc'
      INCLUDE 'interaction_names.inc'

C  Common Block Declarations

      COMMON           /sn_model/
     +                 sn_model_switch
      INTEGER          sn_model_switch

C      COMMON           /sn_interaction_names/
C     +                 interaction_name
C      CHARACTER*57     interaction_name(interaction_total)

C      COMMON           /int_class/
C     +                 list_es, list_cc, list_nc
C      INTEGER          list_es(es_size),
C     +                 list_cc(cc_size),
C     +                 list_nc(nc_size)

C      COMMON           /nu_class/
C     +                 list_nu_e, list_nu_eb,
C     +                 list_nu_mu, list_nu_mub,
C     +                 list_nu_tau, list_nu_taub
C      INTEGER          list_nu_e(nu_e_size),
C     +                 list_nu_eb(nu_eb_size),
C     +                 list_nu_mu(nu_mu_size),
C     +                 list_nu_mub(nu_mub_size),
C     +                 list_nu_tau(nu_tau_size),
C     +                 list_nu_taub(nu_taub_size)

C      COMMON           /reg_class/
C     +                 list_d2o, list_h2o, list_av
C      INTEGER          list_d2o(d2o_size),
C     +                 list_h2o(h2o_size),
C     +                 list_av(av_size)

* output COMMON blocks

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
     +                 interaction_lab(interaction_total),
     +                 interaction_Pb(interaction_total)

      COMMON          /interaction_arrays/
     +                interaction_all,
     +                interaction_es, 
     +                interaction_cc,
     +                interaction_nc,
     +                interaction_nu_e, 
     +                interaction_nu_eb, 
     +                interaction_nu_mu,
     +                interaction_d2o,
     +                interaction_h2o,
     +                interaction_av,
     +                interaction_LAB

C  Global Variable Declarations

      INTEGER          id_interaction(*) ! input
      REAL*8           count_pois(*),    ! input
     +                 count_part(*)     ! input
      INTEGER          icount_pois(*),   ! input
     +                 icount_part(*)    ! input
      INTEGER          exp_id            ! input

C  Local Variable Declarations

      INTEGER          i, id, start
      CHARACTER*57     int_name
      CHARACTER*30     cmodel
      REAL*8           sum_all,
     +                 sum_nu_e, sum_nu_eb, sum_nu_mu,
     +                 sum_cc, sum_nc, sum_es,
     +                 sum_d2o, sum_h2o, sum_av, sum_LAB
      INTEGER          isum_all,
     +                 isum_nu_e, isum_nu_eb, isum_nu_mu,
     +                 isum_cc, isum_nc, isum_es,
     +                 isum_d2o, isum_h2o, isum_av, isum_LAB

***************************************************************************

**************************************************************************
*  Initializations
**************************************************************************

* initialize interaction arrays

      DO i = 1, interaction_total
         interaction_all(i) = 0

         interaction_es(i) = 0
         interaction_cc(i) = 0
         interaction_nc(i) = 0

         interaction_nu_e(i) = 0
         interaction_nu_eb(i) = 0
         interaction_nu_mu(i) = 0

         interaction_d2o(i) = 0
         interaction_h2o(i) = 0
         interaction_av(i) = 0
         interaction_LAB(i) = 0
         interaction_Pb(i) = 0

      END DO

* fill interaction arrays with interaction id values

*******
* SNO
*******

      IF (exp_id .EQ. 1) THEN

       DO i = 1, interaction_total
          interaction_all(i) = i
       END DO

       DO i = 1, es_size_SNO
          interaction_es(i) = list_es_SNO(i)
       END DO

       DO i = 1, cc_size_SNO
          interaction_cc(i) = list_cc_SNO(i)
       END DO

       DO i = 1, nc_size_SNO
          interaction_nc(i) = list_nc_SNO(i)
       END DO

       DO i = 1, nu_e_size_SNO
          interaction_nu_e(i) = list_nu_e_SNO(i)
       END DO

       DO i = 1, nu_eb_size_SNO
          interaction_nu_eb(i) = list_nu_eb_SNO(i)
       END DO

       DO i = 1, nu_mu_size_SNO
          interaction_nu_mu(i) = list_nu_mu_SNO(i)
       END DO
       start = nu_mu_size_SNO
       DO i = 1, nu_mub_size_SNO
          interaction_nu_mu(start+i) = list_nu_mub_SNO(i)
       END DO
       start = nu_mu_size_SNO + nu_mub_size_SNO
       DO i = 1, nu_tau_size_SNO
          interaction_nu_mu(start+i) = list_nu_tau_SNO(i)
       END DO
       start = nu_mu_size_SNO + nu_mub_size_SNO + nu_tau_size_SNO
       DO i = 1, nu_taub_size_SNO
          interaction_nu_mu(start+i) = list_nu_taub_SNO(i)
       END DO

       DO i = 1, d2o_size
          interaction_d2o(i) = list_d2o(i)
       END DO

       DO i = 1, h2o_size
          interaction_h2o(i) = list_h2o(i)
       END DO

       DO i = 1, av_size
          interaction_av(i) = list_av(i)
       END DO

      ENDIF

********
* SNO+
********

      IF (exp_id .EQ. 2) THEN

       DO i = 1, es_size_SNOP
          interaction_all(i) = list_es_SNOP(i)
       END DO
       start = es_size_SNOP
       DO i = 1, cc_size_SNOP
          interaction_all(i+start) = list_cc_SNOP(i)
       END DO
       start = es_size_SNOP + cc_size_SNOP
       DO i = 1, nc_size_SNOP
          interaction_all(i+start) = list_nc_SNOP(i)
       END DO

       DO i = 1, es_size_SNOP
          interaction_es(i) = list_es_SNOP(i)
       END DO

       DO i = 1, cc_size_SNOP
          interaction_cc(i) = list_cc_SNOP(i)
       END DO

       DO i = 1, nc_size_SNOP
          interaction_nc(i) = list_nc_SNOP(i)
       END DO

       DO i = 1, nu_e_size_SNOP
          interaction_nu_e(i) = list_nu_e_SNOP(i)
       END DO

       DO i = 1, nu_eb_size_SNOP
          interaction_nu_eb(i) = list_nu_eb_SNOP(i)
       END DO

       DO i = 1, nu_mu_size_SNOP
          interaction_nu_mu(i) = list_nu_mu_SNOP(i)
       END DO
       start = nu_mu_size_SNOP
       DO i = 1, nu_mub_size_SNOP
          interaction_nu_mu(start+i) = list_nu_mub_SNOP(i)
       END DO
       start = nu_mu_size_SNOP + nu_mub_size_SNOP
       DO i = 1, nu_tau_size_SNOP
          interaction_nu_mu(start+i) = list_nu_tau_SNOP(i)
       END DO
       start = nu_mu_size_SNOP + nu_mub_size_SNOP + nu_tau_size_SNOP
       DO i = 1, nu_taub_size_SNOP
          interaction_nu_mu(start+i) = list_nu_taub_SNOP(i)
       END DO

       DO i = 1, h2o_size
          interaction_h2o(i) = list_h2o(i)
       END DO

       DO i = 1, av_size
          interaction_av(i) = list_av(i)
       END DO

       DO i = 1, LAB_size
          interaction_LAB(i) = list_LAB(i)
       END DO

      ENDIF

********
* HALO
********

      IF (exp_id .EQ. 3) THEN
       
       DO i = 1, es_size_HALO
          interaction_all(i) = list_es_HALO(i)
          interaction_Pb(i) = list_es_HALO(i)
       END DO
       start = es_size_HALO
       DO i = 1, cc_size_HALO
          interaction_all(i+start) = list_cc_HALO(i)
          interaction_Pb(i+start) = list_cc_HALO(i)
       END DO
       start = es_size_HALO + cc_size_HALO
       DO i = 1, nc_size_HALO
          interaction_all(i+start) = list_nc_HALO(i)
          interaction_Pb(i+start) = list_nc_HALO(i)
       END DO
    
       DO i = 1, es_size_HALO
          interaction_es(i) = list_es_HALO(i)
       END DO

       DO i = 1, cc_size_HALO
          interaction_cc(i) = list_cc_HALO(i)
       END DO

       DO i = 1, nc_size_HALO
          interaction_nc(i) = list_nc_HALO(i)
       END DO

       DO i = 1, nu_e_size_HALO
          interaction_nu_e(i) = list_nu_e_HALO(i)
       END DO

       DO i = 1, nu_eb_size_HALO
          interaction_nu_eb(i) = list_nu_eb_HALO(i)
       END DO

       DO i = 1, nu_mu_size_HALO
          interaction_nu_mu(i) = list_nu_mu_HALO(i)
       END DO
       start = nu_mu_size_HALO
       DO i = 1, nu_mub_size_HALO
          interaction_nu_mu(start+i) = list_nu_mub_HALO(i)
       END DO
       start = nu_mu_size_HALO + nu_mub_size_HALO
       DO i = 1, nu_tau_size_HALO
          interaction_nu_mu(start+i) = list_nu_tau_HALO(i)
       END DO
       start = nu_mu_size_HALO + nu_mub_size_HALO + nu_tau_size_HALO
       DO i = 1, nu_taub_size_HALO
          interaction_nu_mu(start+i) = list_nu_taub_HALO(i)
       END DO

      ENDIF

      WRITE (*,*) 
      WRITE (*,*) 'Writing summary file...'

      OPEN (unit = 10, file = 'summary.out', status = 'unknown')

***************************************************************************
* heading
***************************************************************************

      WRITE (10,1)
1     FORMAT ('                      Supernova Neutrino Summary')
      WRITE (10,2)
2     FORMAT ('                      ==============================')

      WRITE (10,*)

      IF (sn_model_switch .EQ. 1) THEN
         cmodel = 'Burrows et al. (1992)'
      ELSEIF (sn_model_switch .EQ. 2) THEN
         cmodel = 'Beacom & Vogel (1998)'
      ELSEIF (sn_model_switch .EQ. 3) THEN
         cmodel = 'Bruenn15'
      ELSEIF (sn_model_switch .EQ. 4) THEN
         cmodel = 'Bruenn25'
      ELSEIF (sn_model_switch .EQ. 5) THEN
         cmodel = 'Mezzacappa'
      ENDIF

      WRITE (10,5) cmodel
5     FORMAT ('Supernova Model: ',a)

      WRITE (10,*)
**************************************************************************
* list reactions with individual count sums
**************************************************************************
* Unfortunately there is some complexity in picking the correct interaction from 
* interaction_name: if unsure look at interaction_names.inc to see why this is
* all necessary

      DO i = 1, interaction_total
         id = id_interaction(i)
         IF (count_pois(id) .GT. 0.d0) THEN   ! only write out chosen reactions

* Determine correct interaction name

           IF (id .LT. first_SNOP) THEN
* This takes care of SNO interactions and SNOP water/AV interactions that were 
* present in SNO
              int_name = interaction_name(id)
           ELSEIF (exp_id .EQ. 2) THEN
* This takes care of new interactions in SNOP and interactions in HALO
              int_name = interaction_name(interaction_total_SNO 
     +     + id - first_SNOP + interaction_total_SNOP - LAB_size + 1)
           ELSEIF (exp_id .EQ. 3) THEN
              int_name = interaction_name(interaction_total_SNO
     +     + interaction_total_SNOP + id - first_HALO + 1)
           ENDIF
         WRITE (10,10) int_name, count_pois(id)
10       FORMAT (a,1x,f8.3)
         ENDIF
      END DO

      WRITE (10,*)

**************************************************************************
* sum neutrino counts of various classes of interactions
**************************************************************************

* all interactions

      CALL group_sum1d(interaction_all,count_pois,sum_all,exp_id)
      CALL group_isum1d(interaction_all,icount_pois,isum_all)

      WRITE (10,20) sum_all, isum_all
20    FORMAT ('TOTAL NEUTRINO COUNTS (all reactions)....................
     +', f8.3,' [',i4,']')

      WRITE (10,*)

* sum interactions by neutrino type

      CALL group_sum1d(interaction_nu_e,count_pois,sum_nu_e,exp_id)
      CALL group_sum1d(interaction_nu_eb,count_pois,sum_nu_eb,exp_id)
      CALL group_sum1d(interaction_nu_mu,count_pois,sum_nu_mu,exp_id)

      CALL group_isum1d(interaction_nu_e,icount_pois,isum_nu_e)
      CALL group_isum1d(interaction_nu_eb,icount_pois,isum_nu_eb)
      CALL group_isum1d(interaction_nu_mu,icount_pois,isum_nu_mu)


      IF (sum_all .GT. 0) THEN
       WRITE (10,30) sum_nu_e, isum_nu_e, sum_nu_e / sum_all * 100.d0
       ELSE
       WRITE (10,30) sum_nu_e, isum_nu_e, 0
      ENDIF
30    FORMAT (2x,'nu_e total neutrino counts (all reactions).............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,31) sum_nu_eb, isum_nu_eb, 
     +              sum_nu_eb / sum_all * 100.d0
       ELSE
       WRITE (10,31) sum_nu_eb, isum_nu_eb, 0
      ENDIF  
31    FORMAT (2x,'nu_ebar total neutrino counts (all reactions)..........
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,32) sum_nu_mu, isum_nu_mu, sum_nu_mu / sum_all* 100.d0
       ELSE
       WRITE (10,32) sum_nu_mu, isum_nu_mu, 0
      ENDIF
32    FORMAT (2x,'nu_"mu" total neutrino counts (all reactions)..........
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      WRITE (10,*)

* sum interactions by interaction type

      CALL group_sum1d(interaction_cc,count_pois,sum_cc,exp_id)
      CALL group_sum1d(interaction_nc,count_pois,sum_nc,exp_id)
      CALL group_sum1d(interaction_es,count_pois,sum_es,exp_id)

      CALL group_isum1d(interaction_cc,icount_pois,isum_cc)
      CALL group_isum1d(interaction_nc,icount_pois,isum_nc)
      CALL group_isum1d(interaction_es,icount_pois,isum_es)

      IF (sum_all .GT. 0) THEN
       WRITE (10,40) sum_cc, isum_cc, sum_cc / sum_all * 100.d0
      ELSE
       WRITE (10,40) sum_cc, isum_cc, 0
      ENDIF
40    FORMAT (2x,'CC total neutrino counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,41) sum_nc, isum_nc, sum_nc / sum_all * 100.d0
      ELSE
       WRITE (10,41) sum_nc, isum_nc, 0
      ENDIF       
41    FORMAT (2x,'NC total neutrino counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,42) sum_es, isum_es, sum_es / sum_all * 100.d0
      ELSE
       WRITE (10,42) sum_es, isum_es, 0
      ENDIF
42    FORMAT (2x,'ES total neutrino counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      WRITE (10,*)

* sum interactions by region type

      IF (exp_id .NE. 3) THEN

       CALL group_sum1d(interaction_h2o,count_pois,sum_h2o,exp_id)
       CALL group_sum1d(interaction_av,count_pois,sum_av,exp_id)
 
       CALL group_isum1d(interaction_h2o,icount_pois,isum_h2o)
       CALL group_isum1d(interaction_av,icount_pois,isum_av)
 
       IF (exp_id .EQ. 1) THEN

        CALL group_sum1d(interaction_d2o,count_pois,sum_d2o,exp_id)
        CALL group_isum1d(interaction_d2o,icount_pois,isum_d2o)

        ELSEIF (exp_id .EQ. 2) THEN

        CALL group_sum1d(interaction_LAB,count_pois,sum_LAB,exp_id)
        CALL group_isum1d(interaction_LAB,icount_pois,isum_LAB)
  
       ENDIF

       IF (sum_all .GT. 0) THEN
        WRITE (10,50) sum_h2o, isum_h2o, sum_h2o / sum_all * 100.d0
       ELSE
        WRITE (10,50) sum_h2o, isum_h2o, 0
       ENDIF
50     FORMAT (2x,'H2O total neutrino counts (all reactions)..............
     +.. ', f8.3,' [',i4,'] (', f5.1,'%)')


       IF (exp_id .EQ. 1) THEN
        IF (sum_all .GT. 0) THEN
         WRITE (10,51) sum_d2o, isum_d2o, sum_d2o / sum_all * 100.d0
        ELSE
         WRITE (10,51) sum_d2o, isum_d2o, 0
        ENDIF
51      FORMAT (2x,'D2O total neutrino counts (all reactions)..............
     +... ', f8.3,' [',i4,'] (', f5.1,'%)')
       ENDIF

       IF (sum_all .GT. 0) THEN
        WRITE (10,52) sum_av, isum_av, sum_av / sum_all * 100.d0
       ELSE
        WRITE (10,52) sum_av, isum_av, 0
       ENDIF
52     FORMAT (2x,'AV total neutrino counts (all reactions)...............
     +.. ', f8.3,' [',i4,'] (', f5.1,'%)')

       IF (exp_id .EQ. 2) THEN
        IF (sum_all .GT. 0) THEN
         WRITE (10,53) sum_LAB, isum_LAB, sum_LAB / sum_all * 100.d0
        ELSE
         WRITE (10,53) sum_LAB, isum_LAB, 0
        ENDIF
53      FORMAT (2x,'LAB total neutrino counts (all reactions)...............
     +... ', f8.3,' [',i4,'] (', f5.1,'%)')
       ENDIF

      ENDIF

**************************************************************************
* sum particle counts of various classes of interactions
* (use same variables for ease)
**************************************************************************

      WRITE (10,*)

* all interactions

      CALL group_sum1d(interaction_all,count_part,sum_all,exp_id)
      CALL group_isum1d(interaction_all,icount_part,isum_all)

      WRITE (10,100) sum_all, isum_all
100   FORMAT ('TOTAL PARTICLE COUNTS (all reactions)....................
     +', f8.3,' [',i4,']')

      WRITE (10,*)

* sum interactions by neutrino type

      CALL group_sum1d(interaction_nu_e,count_part,sum_nu_e,exp_id)
      CALL group_sum1d(interaction_nu_eb,count_part,sum_nu_eb,exp_id)
      CALL group_sum1d(interaction_nu_mu,count_part,sum_nu_mu,exp_id)

      CALL group_isum1d(interaction_nu_e,icount_part,isum_nu_e)
      CALL group_isum1d(interaction_nu_eb,icount_part,isum_nu_eb)
      CALL group_isum1d(interaction_nu_mu,icount_part,isum_nu_mu)

      IF (sum_all .GT. 0) THEN
        WRITE (10,110) sum_nu_e, isum_nu_e, sum_nu_e / sum_all * 100.d0
      ELSE
       WRITE (10,110) sum_nu_e, isum_nu_e, 0
      ENDIF
110   FORMAT (2x,'nu_e total particle counts (all reactions).............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,111) sum_nu_eb, isum_nu_eb, sum_nu_eb / sum_all* 100.d0
      ELSE
       WRITE (10,111) sum_nu_eb, isum_nu_eb, 0
      ENDIF
111   FORMAT (2x,'nu_ebar total particle counts (all reactions)..........
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
      WRITE (10,112) sum_nu_mu, isum_nu_mu, sum_nu_mu / sum_all * 100.d0
      ELSE
       WRITE (10,112) sum_nu_mu, isum_nu_mu, 0
      ENDIF
112   FORMAT (2x,'nu_"mu" total particle counts (all reactions)..........
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      WRITE (10,*)

* sum interactions by interaction type

      CALL group_sum1d(interaction_cc,count_part,sum_cc,exp_id)
      CALL group_sum1d(interaction_nc,count_part,sum_nc,exp_id)
      CALL group_sum1d(interaction_es,count_part,sum_es,exp_id)

      CALL group_isum1d(interaction_cc,icount_part,isum_cc)
      CALL group_isum1d(interaction_nc,icount_part,isum_nc)
      CALL group_isum1d(interaction_es,icount_part,isum_es)

      IF (sum_all .GT. 0) THEN
       WRITE (10,120) sum_cc, isum_cc, sum_cc / sum_all * 100.d0
      ELSE
       WRITE (10,120) sum_cc, isum_cc, 0
      ENDIF
120   FORMAT (2x,'CC total particle counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,121) sum_nc, isum_nc, sum_nc / sum_all * 100.d0
      ELSE
       WRITE (10,121) sum_nc, isum_nc, 0
      ENDIF
121   FORMAT (2x,'NC total particle counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      IF (sum_all .GT. 0) THEN
       WRITE (10,122) sum_es, isum_es, sum_es / sum_all * 100.d0
      ELSE
       WRITE (10,122) sum_es, isum_es, 0
      ENDIF
122   FORMAT (2x,'ES total particle counts (all reactions)...............
     +. ', f8.3,' [',i4,'] (', f5.1,'%)')

      WRITE (10,*)

* sum interactions by region type

      IF (exp_id .NE. 3) THEN

       CALL group_sum1d(interaction_h2o,count_part,sum_h2o,exp_id)
       CALL group_sum1d(interaction_av,count_part,sum_av,exp_id)

       CALL group_isum1d(interaction_h2o,icount_part,isum_h2o)
       CALL group_isum1d(interaction_av,icount_part,isum_av)

       IF (exp_id .EQ. 1) THEN

       CALL group_sum1d(interaction_d2o,count_part,sum_d2o,exp_id)
       CALL group_isum1d(interaction_d2o,icount_part,isum_d2o)

       ELSEIF (exp_id .EQ. 2) THEN

       CALL group_sum1d(interaction_LAB,count_part,sum_LAB,exp_id)
       CALL group_isum1d(interaction_LAB,icount_part,isum_LAB)

       ENDIF

       IF (sum_all .GT. 0) THEN
        WRITE (10,130) sum_h2o, isum_h2o, sum_h2o / sum_all * 100.d0
       ELSE
        WRITE (10,130) sum_h2o, isum_h2o, 0
       ENDIF
130    FORMAT (2x,'H2O total particle counts (all reactions)..............
     +.. ', f8.3,' [',i4,'] (', f5.1,'%)')

       IF (exp_id .EQ. 1) THEN
        IF (sum_all .GT. 0) THEN 
         WRITE (10,131) sum_d2o, isum_d2o, sum_d2o / sum_all * 100.d0
        ELSE
         WRITE (10,131) sum_d2o, isum_d2o, sum_d2o / sum_all * 100.d0
        ENDIF
131     FORMAT (2x,'D2O total particle counts (all reactions)..............
     +... ', f8.3,' [',i4,'] (', f5.1,'%)')
       ENDIF

       IF (sum_all .GT. 0) THEN
        WRITE (10,132) sum_av, isum_av, sum_av / sum_all * 100.d0
       ELSE
        WRITE (10,132) sum_av, isum_av, 0
       ENDIF
132    FORMAT (2x,'AV total particle counts (all reactions)...............
     +.. ', f8.3,' [',i4,'] (', f5.1,'%)')

       IF (exp_id .EQ. 2) THEN
        IF (sum_all .GT. 0) THEN
         WRITE (10,133) sum_LAB, isum_LAB, sum_LAB / sum_all * 100.d0
        ELSE
         WRITE (10,133) sum_LAB, isum_LAB, 0
        ENDIF
133     FORMAT (2x,'LAB total particle counts (all reactions)................
     +... ', f8.3,' [',i4,'] (', f5.1,'%)')
       ENDIF

      ENDIF

       CLOSE (unit = 10)

      RETURN
      END 
