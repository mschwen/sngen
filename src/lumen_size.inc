****************************************************************************
*  Filename: lumen_size.inc
****************************************************************************

      INTEGER          data_size
      PARAMETER        (data_size = 1000)  ! enough to fit any of these

      INTEGER          lum_Bruenn15_entries, lum_Bruenn25_entries,
     +                 lum_Mezz13_entries
      PARAMETER        (lum_Bruenn15_entries = 460)  ! Bruenn 15 Ms
      PARAMETER        (lum_Bruenn25_entries = 999)  ! Bruenn 25 Ms
      PARAMETER        (lum_Mezz13_entries = 665)    ! Mezzacappa 13 Ms

      INTEGER          Eavg_Bruenn15_entries, Eavg_Bruenn25_entries,
     +                 Eavg_Mezz13_entries
      PARAMETER        (Eavg_Bruenn15_entries = 460) ! Bruenn 15 Ms
      PARAMETER        (Eavg_Bruenn25_entries = 999) ! Bruenn 25 Ms
      PARAMETER        (Eavg_Mezz13_entries = 675)   ! Mezzacappa 13 Ms

