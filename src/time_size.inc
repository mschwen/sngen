****************************************************************************
*  Filename: time_size.inc
*    o time_size   = (tstop - tstart) / tstep + 1
****************************************************************************

      INTEGER time_size, time_d_size

* burst: [(1000 - (-100)) / 1] + 1 = 1101
* cooling: [(50000 - 1010) / 10] = 4890  ! default
* cooling: [(50000 - 1020) / 20] = 2449  ! new Mac OS X setting
* --> total = 1101 + 4890 = 5991  ! default
* --> total = 1101 + 2449 = 3550  ! new Mac OS X setting
c      PARAMETER (time_size = 6001)  ! enough for Burrows
      PARAMETER (time_size = 4019)  ! enough for BeacomVogel
c      PARAMETER (time_size = 440)   ! Bruenn 15Ms
c      PARAMETER (time_size = 190)   ! Bruenn 25 Ms
c      PARAMETER (time_size = 665)   ! Mezzacappa 13 Ms


* enough space to analyse m_nu = 100 eV

c      PARAMETER (time_d_size = time_size + 0)
c      PARAMETER (time_d_size = time_size + 400000)
      PARAMETER (time_d_size = time_size + 20000)

