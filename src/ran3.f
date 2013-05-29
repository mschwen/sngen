************************************************************************
*  Filename: ran3.f
************************************************************************
*
*  Random number generator from Numerical Recipes (no modifications)
*  "ran3" takes a -ve number to initialize a set of random
*  numbers (only needed once). After that, the argument that gets
*  passed in is "1" and does not initialize the function again.
*
************************************************************************

      DOUBLE PRECISION FUNCTION ran3(idum)

      IMPLICIT none

C   Global Variable Declarations

      INTEGER idum

C   Local Variable Declarations 

      INTEGER MBIG,MSEED,MZ
      DOUBLE PRECISION FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)

      SAVE iff,inext,inextp,ma
      DATA iff /0/

*************************************************************************

      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj

        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
c        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)

      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3 = mj * FAC       

      IF (ran3 .LT. 0. .OR. ran3 .GT. 1.) THEN
         WRITE (*,*) '*** ran3.f: something went wrong!'
         WRITE (*,*) '*** ran3.f: ran3 =', ran3
      ENDIF

      RETURN
      END
