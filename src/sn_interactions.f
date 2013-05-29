****************************************************************************
*  Filename: sn_interactions.f
****************************************************************************
*
*  Covers the list of all possible supernova interactions based on experiment 
*  (56 for SNO, 52 for SNO+, 8 for HALO); additional choices
*  encompass classes of interactions (ES, CC, NC) as well as
*  the various neutrino types.  Of course, there's also an option 
*  to choose all reactions.
*
****************************************************************************

      SUBROUTINE sn_interactions(interaction_id,exp_id)

      IMPLICIT none

      INCLUDE 'interaction_total.inc'
      INCLUDE 'interaction_names.inc'
      INCLUDE 'sn_classes.inc'

C  Common Block Declarations

      INTEGER pick_group_es, pick_group_cc, pick_group_nc,
     +        pick_group_nu_e, pick_group_nu_eb, pick_group_nu_mu,
     +        pick_group_d2o, pick_group_h2o, pick_group_LAB,
     +        pick_group_av
      COMMON /interaction_group/ 
     + pick_group_es, pick_group_cc, pick_group_nc,
     + pick_group_nu_e, pick_group_nu_eb, pick_group_nu_mu,
     + pick_group_d2o, pick_group_h2o, pick_group_LAB,
     + pick_group_av

C  Local Variable Declarations

      INTEGER stop
      INTEGER interaction
      INTEGER index_start
      LOGICAL done
      LOGICAL check
      INTEGER i, j, k
      INTEGER exp_id

****************************************************************************

****************************************************************************
*  Initializations
****************************************************************************

C  Global Variable Declarations (out)

      INTEGER interaction_id(interaction_total)

****************************************************************************
 
C       WRITE(*,*) "Choose experiment (SNO = 1, SNOP = 2, HALO = 3):"
C       READ(*,*) exp_id

C       IF (exp_id .EQ. 1) THEN
C        interaction_total = interaction_total_SNO
C       ENDIF
C       IF (exp_id .EQ. 2) THEN
C         interaction_total = interaction_total_SNOP
C       ENDIF
C       IF (exp_id .EQ. 3) THEN
C         interaction_total = interaction_total_HALO
C       ENDIF

      DO 2 i = 1, interaction_total
         interaction_id(i) = 0
2     CONTINUE

      interaction = 0

      stop = 111

      pick_group_es    = 0
      pick_group_cc    = 0
      pick_group_nc    = 0
      pick_group_nu_e  = 0
      pick_group_nu_eb = 0
      pick_group_nu_mu = 0
      pick_group_d2o   = 0
      pick_group_h2o   = 0
      pick_group_av    = 0
      pick_group_LAB   = 0

****************************************************************************
c    List interactions based on experiment

      WRITE (*,*) "The interactions for SN neutrinos are as follows:"
      WRITE (*,*) 

      IF (exp_id .EQ. 1) THEN
      i = 0
      DO 5 i = 1, interaction_total_SNO
      WRITE(*,*) interaction_name(i)
5     CONTINUE

      ELSEIF (exp_id .EQ. 2) THEN
      i = 0
      DO 6 i = 1, interaction_total_SNOP
         WRITE(*,*) interaction_name(i + first_SNOP - 1)
6     CONTINUE
	
      ELSEIF (exp_id .EQ. 3) THEN
      i = 0
      DO 7 i = 1, interaction_total_HALO
         WRITE(*,*) interaction_name(i + interaction_total_SNO
     +                                 + interaction_total_SNOP)
7     CONTINUE
      
      ENDIF 
   

       WRITE (*,*)
       PRINT *, "100) Choose all reactions and Proceed"
       WRITE (*,*)

       PRINT *, "101) Choose all ES reactions and Proceed"

       PRINT *, "102) Choose all CC reactions and Proceed"
       PRINT *, "103) Choose all NC reactions and Proceed"
       WRITE (*,*)
       PRINT *, "104) Choose all nu_e reactions and Proceed"
       PRINT *, "105) Choose all nu_e_bar reactions and Proceed"
       PRINT *, '106) Choose all nu_"mu" reactions and Proceed'
       WRITE (*,*)

        IF (exp_id .EQ. 1) THEN
         PRINT *, '107) Choose all D2O reactions and Proceed'
         ELSEIF (exp_id .EQ. 2) THEN
         PRINT *, '108) Choose all LAB reactions and Proceed'
        ENDIF 
        IF (exp_id .NE. 3) THEN
         PRINT *, '109) Choose all H2O reactions and Proceed'
         PRINT *, '110) Choose all AV  reactions and Proceed'
         WRITE (*,*)
        ENDIF
       PRINT *, "111) Proceed"
       WRITE (*,*)

       IF (exp_id .EQ. 1) THEN
        PRINT *, "Please select the interaction(s) [1-56]"
       ELSEIF (exp_id .EQ. 2) THEN
        PRINT *, "Please select the interaction(s) [57-70]"
       ELSEIF (exp_id .EQ. 3) THEN
        PRINT *, "Please select the interaction(s) [75-96]"
       ENDIF

       PRINT *, "or groups of interactions [100-110]"
       PRINT *, "Choose [111] to proceed..."
  

************************************************************************
C  SNO
************************************************************************
      IF (exp_id .EQ. 1) THEN

      done = .false.
      j = 0

      DO WHILE (.NOT. done)

         j = j + 1
         k = 0

         READ (*,*) interaction

* neutrino type reactions
         DO 10 k = 1, nu_e_size_SNO
            IF (interaction .EQ. list_nu_e_SNO(k)) pick_group_nu_e = 1
10       CONTINUE

         DO 11 k = 1, nu_eb_size_SNO
            IF (interaction .EQ. list_nu_eb_SNO(k)) THEN
            pick_group_nu_eb = 1
            ENDIF
11       CONTINUE

         DO 12 k = 1, nu_mu_size_SNO
            IF (interaction .EQ. list_nu_mu_SNO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
12       CONTINUE

         DO 13 k = 1, nu_mub_size_SNO
            IF (interaction .EQ. list_nu_mub_SNO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
13       CONTINUE

         DO 14 k = 1, nu_tau_size_SNO
            IF (interaction .EQ. list_nu_tau_SNO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
14       CONTINUE

         DO 15 k = 1, nu_taub_size_SNO
            IF (interaction .EQ. list_nu_taub_SNO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
15       CONTINUE

* region type reactions
         DO 20 k = 1, d2o_size
            IF (interaction .EQ. list_d2o(k)) pick_group_d2o = 1
20       CONTINUE

         DO 21 k = 1, h2o_size
            IF (interaction .EQ. list_h2o(k)) pick_group_h2o = 1
21       CONTINUE

         DO 22 k = 1, av_size
            IF (interaction .EQ. list_av(k)) pick_group_av = 1
22       CONTINUE

* interaction type reactions
         DO 30 k = 1, cc_size_SNO
            IF (interaction .EQ. list_cc_SNO(k)) pick_group_cc = 1
30       CONTINUE

         DO 31 k = 1, nc_size_SNO
            IF (interaction .EQ. list_nc_SNO(k)) pick_group_nc = 1
31       CONTINUE

         DO 32 k = 1, es_size_SNO
            IF (interaction .EQ. list_es_SNO(k)) pick_group_es = 1
32       CONTINUE

* have all the desired interactions
         IF (interaction .EQ. stop) THEN
            done = .true.
         ELSE
          IF (interaction .LT. 1 .OR. interaction .GT. stop .OR.
     +       (interaction .LT. 99 .AND. interaction .GT. last_SNO)) THEN
               PRINT *, "Try again"
               done = .false.

************************************************************************
* choosing interactions one by one
************************************************************************
 
           ELSEIF (interaction .LE. last_SNO) THEN
               interaction_id(j) = interaction
               done = .false.
            
************************************************************************
* choose all interactions
************************************************************************

            ELSEIF (interaction .EQ. 100) THEN
               DO 40 k = 1, interaction_total
                  interaction_id(k) = k
40             CONTINUE
               done = .true.

************************************************************************
* choose reactions by interaction type
************************************************************************

* all ES reactions

            ELSEIF (interaction .EQ. 101) THEN
               DO 45 k = 1, es_size_SNO
                  interaction_id(k) = list_es_SNO(k)
45             CONTINUE
               done = .true.

* all CC reactions

            ELSEIF (interaction .EQ. 102) THEN
               DO 46 k = 1, cc_size_SNO
                  interaction_id(k) = list_cc_SNO(k)
46             CONTINUE
               done = .true.

* all NC reactions

            ELSEIF (interaction .EQ. 103) THEN
               DO 47 k = 1, nc_size_SNO
                  interaction_id(k) = list_nc_SNO(k)
47             CONTINUE
               done = .true.

************************************************************************
* choose reactions by neutrino type
************************************************************************

* all nu_e reactions

            ELSEIF (interaction .EQ. 104) THEN
               DO 50 k = 1, nu_e_size_SNO
                  interaction_id(k) = list_nu_e_SNO(k)
50             CONTINUE
               done = .true.

* all nu_e_bar reactions

            ELSEIF (interaction .EQ. 105) THEN
               DO 60 k = 1, nu_eb_size_SNO
                  interaction_id(k) = list_nu_eb_SNO(k)
60             CONTINUE
               done = .true.

* all nu_"mu" reactions

            ELSEIF (interaction .EQ. 106) THEN
               DO 70 k = 1, nu_mu_size_SNO
                  index_start = 0
                  interaction_id(k+index_start) = list_nu_mu_SNO(k)
70             CONTINUE
               DO 80 k = 1, nu_mub_size_SNO
                  index_start = nu_mu_size_SNO
                  interaction_id(k+index_start) = list_nu_mub_SNO(k)
80             CONTINUE
               DO 90 k = 1, nu_tau_size_SNO
                  index_start = nu_mu_size_SNO + nu_mub_size_SNO
                  interaction_id(k+index_start) = list_nu_tau_SNO(k)
90             CONTINUE
               DO 100 k = 1, nu_taub_size_SNO
                  index_start = nu_mu_size_SNO + nu_mub_size_SNOP 
     +            + nu_tau_size_SNO
                  interaction_id(k+index_start) = list_nu_taub_SNO(k)
100             CONTINUE
               done = .true.

***************************************************************************
* choose reaction by region type
***************************************************************************

* all D2O reactions

            ELSEIF (interaction .EQ. 107) THEN
               DO 110 k = 1, d2o_size
                  interaction_id(k) = list_d2o(k)
110            CONTINUE
               done = .true.

* all H2O reactions

            ELSEIF (interaction .EQ. 109) THEN
               DO 120 k = 1, h2o_size
                  interaction_id(k) = list_h2o(k)
120            CONTINUE
               done = .true.

* all AV reactions

            ELSEIF (interaction .EQ. 110) THEN
               DO 130 k = 1, av_size
                  interaction_id(k) = list_av(k)
130             CONTINUE
               done = .true.

            ENDIF
         ENDIF    ! check if interaction equals "stop"
      END DO      ! do while not done loop

      ENDIF



*****************************************************************
C   SNO+
*****************************************************************


      IF (exp_id .EQ. 2) THEN

      done = .false.
      j = 0

      DO WHILE (.NOT. done)

         j = j + 1
         k = 0

         READ (*,*) interaction

* neutrino type reactions
         DO 131 k = 1, nu_e_size_SNOP
            IF (interaction .EQ. list_nu_e_SNOP(k)) THEN 
            pick_group_nu_e = 1
            ENDIF
131       CONTINUE

         DO 132 k = 1, nu_eb_size_SNOP
            IF (interaction .EQ. list_nu_eb_SNOP(k)) THEN 
            pick_group_nu_eb = 1
            ENDIF
132       CONTINUE

         DO 133 k = 1, nu_mu_size_SNOP
            IF (interaction .EQ. list_nu_mu_SNOP(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
133       CONTINUE

         DO 134 k = 1, nu_mub_size_SNOP
            IF (interaction .EQ. list_nu_mub_SNOP(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
134       CONTINUE

         DO 135 k = 1, nu_tau_size_SNOP
            IF (interaction .EQ. list_nu_tau_SNOP(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
135       CONTINUE

         DO 136 k = 1, nu_taub_size_SNOP
            IF (interaction .EQ. list_nu_taub_SNOP(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
136       CONTINUE

         DO 137 k =1, LAB_size
            IF (interaction .EQ. list_LAB(k)) THEN
            pick_group_LAB = 1
            ENDIF
137       CONTINUE

         DO 138 k = 1, h2o_size
            IF (interaction .EQ. list_h2o(k)) pick_group_h2o = 1
138       CONTINUE

         DO 139 k = 1, av_size
            IF (interaction .EQ. list_av(k)) pick_group_av = 1
139       CONTINUE

* interaction type reactions
         DO 140 k = 1, cc_size_SNOP
            IF (interaction .EQ. list_cc_SNOP(k)) pick_group_cc = 1
140       CONTINUE

         DO 141 k = 1, nc_size_SNOP
            IF (interaction .EQ. list_nc_SNOP(k)) pick_group_nc = 1
141       CONTINUE

         DO 143 k = 1, es_size_SNOP
            IF (interaction .EQ. list_es_SNOP(k)) pick_group_es = 1
143       CONTINUE

* have all the desired interactions
         IF (interaction .EQ. stop) THEN
            done = .true.
         ELSE
            IF (interaction .LT. 1 .OR. interaction .GT. stop .OR.
     +         (interaction .GT. last_SNOP .AND. 
     +          interaction .LT. 100)) THEN
               PRINT *, "Try again"
               done = .false.

************************************************************************
* choosing interactions one by one
************************************************************************

           ELSEIF (interaction .LE. last_SNOP+1) THEN
               check = .true.
               DO 900 k = 1, d2o_size
               IF (interaction .EQ. list_d2o(k)) check = .false.
900            CONTINUE

               IF (check .EQV. .true.) THEN
                interaction_id(j) = interaction
                done = .false.
                ELSE
                PRINT *, "Try again"
                done = .false.
               ENDIF

************************************************************************
* choose all interactions
************************************************************************

            ELSEIF (interaction .EQ. 100) THEN
               index_start = 0
               DO 144 k = 1, es_size_SNOP
                interaction_id(k) = list_es_SNOP(k)
144            CONTINUE
               index_start = es_size_SNOP
               DO 600 k = 1, cc_size_SNOP
                interaction_id(k+index_start) = list_cc_SNOP(k)
600            CONTINUE
               index_start = es_size_SNOP + cc_size_SNOP
               DO 700 k = 1, nc_size_SNOP
                interaction_id(k+index_start) = list_nc_SNOP(k)
700            CONTINUE
               done = .true.

************************************************************************
* choose reactions by interaction type
************************************************************************

* all ES reactions

            ELSEIF (interaction .EQ. 101) THEN
               DO 145 k = 1, es_size_SNOP
                  interaction_id(k) = list_es_SNOP(k)
145             CONTINUE
               done = .true.

* all CC reactions

            ELSEIF (interaction .EQ. 102) THEN
               DO 146 k = 1, cc_size_SNOP
                  interaction_id(k) = list_cc_SNOP(k)
146             CONTINUE
               done = .true.

* all NC reactions

            ELSEIF (interaction .EQ. 103) THEN
               DO 147 k = 1, nc_size_SNOP
                  interaction_id(k) = list_nc_SNOP(k)
147             CONTINUE
               done = .true.

************************************************************************
* choose reactions by neutrino type
************************************************************************

* all nu_e reactions

            ELSEIF (interaction .EQ. 104) THEN
               DO 150 k = 1, nu_e_size_SNOP
                  interaction_id(k) = list_nu_e_SNOP(k)
150             CONTINUE
               done = .true.

* all nu_e_bar reactions

            ELSEIF (interaction .EQ. 105) THEN
               DO 160 k = 1, nu_eb_size_SNOP
                  interaction_id(k) = list_nu_eb_SNOP(k)
160             CONTINUE
               done = .true.

* all nu_"mu" reactions

            ELSEIF (interaction .EQ. 106) THEN
               DO 170 k = 1, nu_mu_size_SNOP
                  index_start = 0
                  interaction_id(k+index_start) = list_nu_mu_SNOP(k)
170             CONTINUE
               DO 180 k = 1, nu_mub_size_SNOP
                  index_start = nu_mu_size_SNOP
                  interaction_id(k+index_start) = list_nu_mub_SNOP(k)
180             CONTINUE
               DO 190 k = 1, nu_tau_size_SNOP
                  index_start = nu_mu_size_SNOP + nu_mub_size_SNOP
                  interaction_id(k+index_start) = list_nu_tau_SNOP(k)
190             CONTINUE
               DO 200 k = 1, nu_taub_size_SNOP
                  index_start = nu_mu_size_SNOP 
     +            + nu_mub_size_SNOP + nu_tau_size_SNOP
                  interaction_id(k+index_start) = list_nu_taub_SNOP(k)
200             CONTINUE
               done = .true.

***************************************************************************
* choose reaction by region type
***************************************************************************

* all LAB reactions

            ELSEIF (interaction .EQ. 108) THEN
               DO 210 k = 1, LAB_size
                  interaction_id(k) = list_LAB(k)
210            CONTINUE
               done = .true.

* all H2O reactions

            ELSEIF (interaction .EQ. 109) THEN
               DO 220 k = 1, h2o_size
                  interaction_id(k) = list_h2o(k)
220            CONTINUE
               done = .true.

* all AV reactions

            ELSEIF (interaction .EQ. 110) THEN
               DO 230 k = 1, av_size
                  interaction_id(k) = list_av(k)
230             CONTINUE
               done = .true.

            ENDIF
         ENDIF    ! check if interaction equals "stop"
      END DO      ! do while not done loop

      ENDIF

*****************************************************************
C HALO
*****************************************************************

      IF (exp_id .EQ. 3) THEN

      done = .false.
      j = 0

      DO WHILE (.NOT. done)

         j = j + 1
         k = 0

         READ (*,*) interaction

* neutrino type reactions
         DO 240 k = 1, nu_e_size_HALO
            IF (interaction .EQ. list_nu_e_HALO(k)) THEN
            pick_group_nu_e = 1
            ENDIF
240       CONTINUE

         DO 241 k = 1, nu_eb_size_HALO
            IF (interaction .EQ. list_nu_eb_HALO(k)) THEN
            pick_group_nu_eb = 1
            ENDIF
241       CONTINUE

         DO 242 k = 1, nu_mu_size_HALO
            IF (interaction .EQ. list_nu_mu_HALO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
242       CONTINUE

         DO 243 k = 1, nu_mub_size_HALO
            IF (interaction .EQ. list_nu_mub_HALO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
243       CONTINUE

         DO 244 k = 1, nu_tau_size_HALO
            IF (interaction .EQ. list_nu_tau_HALO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
244       CONTINUE

         DO 245 k = 1, nu_taub_size_HALO
            IF (interaction .EQ. list_nu_taub_HALO(k)) THEN
            pick_group_nu_mu = 1
            ENDIF
245       CONTINUE

* interaction type reactions
         DO 250 k = 1, cc_size_HALO
            IF (interaction .EQ. list_cc_HALO(k)) pick_group_cc = 1
250       CONTINUE

         DO 251 k = 1, nc_size_HALO
            IF (interaction .EQ. list_nc_HALO(k)) pick_group_nc = 1
251       CONTINUE

* have all the desired interactions
         IF (interaction .EQ. stop) THEN
            done = .true.
         ELSE
           IF (interaction .LT. first_HALO 
     +         .OR. interaction .GT. stop .OR.
     +         (interaction .GT. last_HALO .AND. 
     +         interaction .LT. 100)) THEN
               PRINT *, "Try again"
               done = .false.

************************************************************************
* choosing interactions one by one
************************************************************************

           ELSEIF (interaction .LE. last_HALO) THEN
               interaction_id(j) = interaction
               done = .false.

************************************************************************
* choose all interactions
************************************************************************

            ELSEIF (interaction .EQ. 100) THEN
               index_start = 0
               DO 252 k = 1, es_size_HALO
                interaction_id(k+index_start) = list_es_HALO(k)
252            CONTINUE
               index_start = es_size_HALO
               DO 260 k = 1, cc_size_HALO
                interaction_id(k+index_start) = list_cc_HALO(k)
260            CONTINUE
               index_start = es_size_HALO + cc_size_HALO
               DO 800 k = 1, nc_size_HALO
                interaction_id(k+index_start) = list_nc_HALO(k)
800            CONTINUE
               done = .true.

************************************************************************
* choose reactions by interaction type
************************************************************************

* all ES reactions

            ELSEIF (interaction .EQ. 101) THEN
               DO 265 k = 1, es_size_HALO
                  interaction_id(k) = list_es_HALO(k)
265             CONTINUE
               done = .true.

* all CC reactions

            ELSEIF (interaction .EQ. 102) THEN
               DO 270 k = 1, cc_size_HALO
                  interaction_id(k) = list_cc_HALO(k)
270             CONTINUE
               done = .true.

* all NC reactions

            ELSEIF (interaction .EQ. 103) THEN
               DO 280 k = 1, nc_size_HALO
                  interaction_id(k) = list_nc_HALO(k)
280             CONTINUE
               done = .true.

************************************************************************
* choose reactions by neutrino type
************************************************************************

* all nu_e reactions

            ELSEIF (interaction .EQ. 104) THEN
               DO 290 k = 1, nu_e_size_HALO
                  interaction_id(k) = list_nu_e_HALO(k)
290             CONTINUE
               done = .true.

* all nu_e_bar reactions

            ELSEIF (interaction .EQ. 105) THEN
               DO 300 k = 1, nu_eb_size_HALO
                  interaction_id(k) = list_nu_eb_HALO(k)
300             CONTINUE
               done = .true.

* all nu_"mu" reactions

            ELSEIF (interaction .EQ. 106) THEN
               DO 370 k = 1, nu_mu_size_HALO
                  index_start = 0
                  interaction_id(k+index_start) = list_nu_mu_HALO(k)
370             CONTINUE
               DO 380 k = 1, nu_mub_size_HALO
                  index_start = nu_mu_size_HALO
                  interaction_id(k+index_start) = list_nu_mub_HALO(k)
380             CONTINUE
               DO 390 k = 1, nu_tau_size_HALO
                  index_start = nu_mu_size_HALO + nu_mub_size_HALO
                  interaction_id(k+index_start) = list_nu_tau_HALO(k)
390             CONTINUE
               DO 400 k = 1, nu_taub_size_HALO
                  index_start = nu_mu_size_HALO + nu_mub_size_HALO + 
     +            nu_tau_size_HALO
                  interaction_id(k+index_start) = list_nu_taub_HALO(k)
400             CONTINUE
               done = .true.
 
            ENDIF
         ENDIF    ! check if interaction equals "stop"
      END DO      ! do while not done loop

      ENDIF


* determine whether or not all of one nu type has been selected
* nu_type count rate stuff

      IF (interaction .EQ. interaction_total + 100 .OR. 
     +    interaction .EQ. interaction_total + 104) THEN
          pick_group_nu_e = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR. 
     +    interaction .EQ. interaction_total + 105) THEN
          pick_group_nu_eb = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR. 
     +    interaction .EQ. interaction_total + 106) THEN
          pick_group_nu_mu = 1
      ENDIF

*** event class count rate stuff

      IF (interaction .EQ. interaction_total + 100 .OR.   
     +    interaction .EQ. interaction_total + 101) THEN
          pick_group_es = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR.
     +    interaction .EQ. interaction_total + 102) THEN
          pick_group_cc = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR.
     +    interaction .EQ. interaction_total + 103) THEN
          pick_group_nc = 1
      ENDIF 

*** region countrate/count stuff

      IF (interaction .EQ. interaction_total + 100 .OR.
     +    interaction .EQ. interaction_total + 107) THEN
          pick_group_d2o = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR.
     +    interaction .EQ. interaction_total + 108) THEN
          pick_group_h2o = 1
      ENDIF
      IF (interaction .EQ. interaction_total + 100 .OR.
     +    interaction .EQ. interaction_total + 110) THEN
          pick_group_av = 1
      ENDIF

C      DO 440 i = 1, interaction_total
C         WRITE (*,*) 'Interaction =', interaction_id(i)
C440   CONTINUE

      RETURN
      END
