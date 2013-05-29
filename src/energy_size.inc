***************************************************************************
*  Filename: energy_size.inc
***************************************************************************

      INTEGER energy_size
      INTEGER file_size_Br15
      INTEGER file_size_Br25
      INTEGER file_size_Bur
      INTEGER file_size_Mez13
      INTEGER file_size_large

* Burrows et al. (1992)
c      PARAMETER (energy_size = 1000)   ! E_nu to 100 MeV

* Beacom & Vogel (1998)
      PARAMETER (energy_size = 1700)   ! E_nu to 170 MeV

* Bruenn, Mezzacappa (2000)
c      PARAMETER (energy_size = 20)     ! various

* Size of critical energy in collective flavour transitions

      PARAMETER (file_size_Br15 = 40)
      PARAMETER (file_size_Br25 = 55)
      PARAMETER (file_size_Bur = 51)
      PARAMETER (file_size_Mez13 = 26)
      PARAMETER (file_size_large = 55) ! choose largest of the above


