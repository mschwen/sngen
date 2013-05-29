*file member=EPHEM library=snoman language=fortran77 date=18:Aug:1999
      SUBROUTINE EPHEM (IYEAR,IMONTH,TZERO,ICNS,DCONS)
 
*     SOUDAN astronomy routine modified for use in SNOMAN.
*     Contact: R. Van de Water (Penn)
 
*-------------------------------------------------------------------------
*--
*--  Specification
*--  -------------
*--
*--  ALL VALUES IN SECONDS
*--  Constants for a particular year....
*--            TZERO   : GMST at Jan 0 UT=0 for that particular year
*--            AT      : TAI - UTC
*--            TT      : TDT - TAI
*--            DT      : TDT - UT
*--
*--  Some definitions ( see for example B4,L1 astro almanac )...
*--            GMST    : Greenwich Mean Sidereal Time.
*--            UT      : Universal time counted from 0h midnight, its
*--                        unit is the mean solar day.
*--            UTC     : Coordinated universal time, this is the basis of
*--                        the WWVB signal.
*--            TAI     : International atomic time. Differs from UTC by
*--                        an integer number of seconds.
*--            TDT     : terrestial dynamic time, the basis for ephemerides
*--                        for earth bound  observations.
*--
*--  Author  MAT (SOUDAN)
*--  ------
*--  Modified by R.G. Van de Water (Penn): convert to integer format,
*--                                        update for 1993 and 1994
*--  Modified by W. Frati          (Penn): To use past 1993
*---------------------------------------------------------------------------
      IMPLICIT NONE
 
      INTEGER          IYEAR,IMONTH
      INTEGER          JYEAR
      INTEGER          AT,TT,DT
      INTEGER          ICNS(3)
      DOUBLE PRECISION AT_NSEC,TT_NSEC,DT_NSEC
      DOUBLE PRECISION DCONS(3),TZERO
      DOUBLE PRECISION DAY,YEAR
      DOUBLE PRECISION GMST_TZERO
 
C       ASSIGNS APPROPRIATE CONSTANTS FOR TRANSFORMING BETWEEN TIME SYSTEMS
C       ALL CONSTANTS IN SECONDS
C       TZERO   = GMST AT UT = 0 FOR JAN 0 OF THAT YEAR
C
      JYEAR = IYEAR
      IF(IYEAR .LT. 1990) JYEAR = 1990
      YEAR = IYEAR
      DAY = (IMONTH-1)/12. + YEAR
      TZERO     = GMST_TZERO(IYEAR)
      AT_NSEC   = 0.0D0
      TT        = 32
      TT_NSEC   = 1.84D8
C
C     The following calculation of DT is based on the straight line fit
C     from the data in the Astronomical Almanac Pg K9
C     Data from 1995 on is extrapolated
C     Should be updated as time goes by
C
      DT        = INT(9.51997*(DAY-1990.)/10.25  + 56.43566)
      DT_NSEC   = 0.0D0
C
C     The following definition of AT is off by (UTC-UT)
      AT        = DT - TT
      AT_NSEC   = DT_NSEC - TT_NSEC
      ICNS(1)=  AT
      DCONS(1)=  AT_NSEC
      ICNS(2)=  TT
      DCONS(2)=  TT_NSEC
      ICNS(3)=  DT
      DCONS(3)=  DT_NSEC
c     WRITE(6,*) IYEAR,TZERO,DT,AT
      END
*endfile member=EPHEM
