*************************************************************************
*  Filename: convert_units.inc
*************************************************************************

      DOUBLE PRECISION kpc2m, kpc2mev, cm2m, erg2mev, 
     +                 sec2year, ms2s, ms2ns, sec2ms

**** [Particle Data Book]
      PARAMETER (kpc2m = 3.0856775807d+19)

**** [(actually, inverse MeV) NIST]
      PARAMETER (kpc2mev = 2.4887658959d+31)

**** [convert centimeters to meters]
      PARAMETER (cm2m = 1.d-02)

**** [Particle Data Book: 10d-07 / (10d+06 * 1.602d-19)]
      PARAMETER (erg2mev = 6.2415063631d+05)

**** [Particle Data Book: Tropical year = 31556925.2 seconds]
      PARAMETER (sec2year = 3.168876542d-08)

**** [convert milliseconds to seconds]
      PARAMETER (ms2s = 1.d-03)

**** [convert milliseconds to nanoseconds]
      PARAMETER (ms2ns = 1.d+06)

**** [convert seconds to milliseconds]
      PARAMETER (sec2ms = 1.d+03)
