      INCLUDE 'int_class_size.inc'
      INCLUDE 'nu_class_size.inc'
      INCLUDE 'reg_class_size.inc'

      INTEGER list_es_SNO(es_size_SNO), 
     +        list_cc_SNO(cc_size_SNO),
     +        list_nc_SNO(nc_size_SNO),
     +        list_es_SNOP(es_size_SNOP), 
     +        list_cc_SNOP(cc_size_SNOP),
     +        list_nc_SNOP(nc_size_SNOP),
     +        list_es_HALO(es_size_HALO),
     +        list_cc_HALO(cc_size_HALO),
     +        list_nc_HALO(nc_size_HALO)

      INTEGER list_nu_e_SNO(nu_e_size_SNO),
     +        list_nu_eb_SNO(nu_eb_size_SNO),
     +        list_nu_mu_SNO(nu_mu_size_SNO),
     +        list_nu_mub_SNO(nu_mub_size_SNO),
     +        list_nu_tau_SNO(nu_tau_size_SNO),
     +        list_nu_taub_SNO(nu_taub_size_SNO)

      INTEGER list_nu_e_SNOP(nu_e_size_SNOP), 
     +        list_nu_eb_SNOP(nu_eb_size_SNOP),
     +        list_nu_mu_SNOP(nu_mu_size_SNOP),
     +        list_nu_mub_SNOP(nu_mub_size_SNOP),
     +        list_nu_tau_SNOP(nu_tau_size_SNOP),
     +        list_nu_taub_SNOP(nu_taub_size_SNOP)


      INTEGER list_nu_e_HALO(nu_e_size_HALO),
     +        list_nu_eb_HALO(nu_eb_size_HALO),
     +        list_nu_mu_HALO(nu_mu_size_HALO),
     +        list_nu_mub_HALO(nu_mub_size_HALO),
     +        list_nu_tau_HALO(nu_tau_size_HALO),
     +        list_nu_taub_HALO(nu_taub_size_HALO)

      INTEGER list_d2o(d2o_size), list_h2o(h2o_size), list_av(av_size),
     +        list_LAB(LAB_size)

      INTEGER last_SNO, last_SNOP, last_HALO, first_HALO, first_SNOP

C      COMMON  /int_class/
C     +        list_es_SNO, list_cc_SNO, list_nc_SNO,
C     +        list_es_SNOP, list_cc_SNOP, list_nc_SNOP,
C     +        list_cc_HALO, list_nc_HALO

C      COMMON  /nu_class/
C     +        list_nu_e_SNO, list_nu_eb_SNO,
C     +        list_nu_mu_SNO, list_nu_mub_SNO,
C     +        list_nu_tau_SNO, list_nu_taub_SNO,
C     +        list_nu_e_SNOP, list_nu_eb_SNOP,
C     +        list_nu_mu_SNOP, list_nu_mub_SNOP,
C     +        list_nu_tau_SNOP, list_nu_taub_SNOP,
C     +        list_nu_e_HALO, list_nu_eb_HALO,
C     +        list_nu_mu_HALO, list_nu_mub_HALO,
C     +        list_nu_tau_HALO, list_nu_taub_HALO

C      COMMON  /reg_class/
C     +        list_d2o, list_h2o, list_av, list_LAB

      DATA    list_es_SNO /1,2,10,11,21,22,30,31,39,40,48,49/
      DATA    list_cc_SNO /3,4,5,6,7,12,13,14,15,16,17,18/
      DATA    list_nc_SNO /8,9,19,20,
     +                 23,24,25,26,27,28,29,
     +                 32,33,34,35,36,37,38,
     +                 41,42,43,44,45,46,47,
     +                 50,51,52,53,54,55,56/
    
      DATA    list_es_SNOP /1,10,21,30,39,48,57,58,61,62,66,67,69,70,
     +                      72,73,75,76/
      DATA    list_cc_SNOP /4,5,7,12,14,16,17,59,63,64/
      DATA    list_nc_SNOP /9,20,24,26,28,33,35,37,42,44,46,51,53,55,
     +                     60,65,68,71,74,77/

      DATA    list_es_HALO /78,83,88,91,94,97/
      DATA    list_cc_HALO /79,80,84,85/
      DATA    list_nc_HALO /81,82,86,87,89,90,92,93,95,96,98,99/
 

      DATA    list_nu_e_SNO    /1,2,3,4,5,6,7,8,9/
      DATA    list_nu_eb_SNO   /10,11,12,13,14,15,16,17,18,19,20/
      DATA    list_nu_mu_SNO   /21,22,23,24,25,26,27,28,29/
      DATA    list_nu_mub_SNO  /30,31,32,33,34,35,36,37,38/
      DATA    list_nu_tau_SNO  /39,40,41,42,43,44,45,46,47/
      DATA    list_nu_taub_SNO /48,49,50,51,52,53,54,55,56/

      DATA    list_nu_e_SNOP    /1,4,5,7,9,57,58,59,60/
      DATA    list_nu_eb_SNOP   /10,12,13,14,16,17,20,61,62,63,64,65/
      DATA    list_nu_mu_SNOP   /21,24,26,28,66,67,68/
      DATA    list_nu_mub_SNOP  /30,33,35,37,69,70,71/
      DATA    list_nu_tau_SNOP  /39,42,44,46,72,73,74/
      DATA    list_nu_taub_SNOP /48,51,53,55,75,76,77/

      DATA    list_nu_e_HALO    /78,79,80,81,82/
      DATA    list_nu_eb_HALO   /83,84,85,86,87/
      DATA    list_nu_mu_HALO   /88,89,90/
      DATA    list_nu_mub_HALO  /91,92,93/
      DATA    list_nu_tau_HALO  /94,95,96/
      DATA    list_nu_taub_HALO /97,98,99/

      DATA    list_d2o /2,3,6,8,11,13,15,18,19,22,23,25,27,29,31,32,
     +                  34,36,38,40,41,43,45,47,49,50,52,54,56/
      DATA    list_h2o /1,5,10,12,17,21,26,28,30,35,37,39,44,46,48,
     +                  53,55/
      DATA    list_av  /4,7,9,14,16,20,24,33,42,51/
      DATA    list_LAB /57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,
     +                  72,73,74,75,76,77/

      PARAMETER (last_SNO = 56, last_SNOP = 77, first_SNOP = 57,
     + first_HALO = 78, last_HALO = 99)
