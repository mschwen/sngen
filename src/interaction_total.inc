      INTEGER interaction_total_SNO
      INTEGER interaction_total_SNOP
      INTEGER interaction_total_HALO
      INTEGER interaction_total


* Don't change from 100 as in many places the code assumes
* that all interaction numbers are smaller than the interaction
* total - this was true when only SNO was part of the code.
 
      PARAMETER (interaction_total_SNO = 56)
      PARAMETER (interaction_total_SNOP = 49)
      PARAMETER (interaction_total_HALO = 22)
      PARAMETER (interaction_total = 100)
