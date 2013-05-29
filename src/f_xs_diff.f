*FILE MEMBER=F_XS_DIFF LIBRARY=SNOMAN LANGUAGE=FORTRAN77 DATE=5:Jun:2001
      REAL FUNCTION F_XS_DIFF(XSTYPE,ENU,COS_OUT,EOUT)
 
*     Monte Carlo generation - function to return double differential
*	    cross section dsigma/dTdcos(theta) for neutrons produce via
*           neutral current interaction or elecrons produced via
*	    charged current.
 
*     Contact:  J.R. Klein, Penn
*     Verified:
*     Refereed:
 
*     Parameters:-
*     ==========
 
*     ENU          in      Neutrino energy
*     COS_OUT      in	   Cos(theta) for outgoing particle relative to neutrino
*     EOUT	   in	   Energy of outgoing particle
 
*     Common Block Access:-
*     ===================
*     (ignoring internal and environmental access)
 
*     /mt/            in       to obtain RCONS array
*     /titles_ncxs/   in       for links
 
*     Specification:-
*     =============
 
* o   Return double differential for NC or CC cross sections
*
 
 
*     Assumptions:-
*     ===========
*
*	  The double differential cross sections used here are those produced
*	by Kubodera and Nozawa.  The numbers come from a personal communication
*	from Kubodera to Richard Van de Water.
*         The cross sections range from ENU=2.4 MeV to ENU=170.0 MeV.
*         NSGK cross sections range from ENU=2.0 MeV to ENU=170.0 MeV for CC and
*       ENU=2.3 MeV to  ENU=50 MeV.
 
*     Program Notes:-
*     =============
 
*
*   Cross sections are stored in RCONS, and pointers for relevant indices
*  are also from RCONS.  See descriptions of format in titles_ccxs.doc,
*  titles_ccxb.doc, titles_ncxs.doc and titles_ncxb.doc.
*   Kubodera's cross sections were normalized by dividing by the appropriate
*  total cross section at each neutrino energy.  This is undone here and
*  differential cross section is multiplied by the total cross section before
*  the final differential cross section is calculated, returning a number
*  with the units 10^-44*cm^2/MeV
*   For neutrino energies which are above threshold but below the lowest
*  energy given in Kubodera's table, an extrapolation toward threshold is
*  done using the next two highest energy points in the table.  A parabola
*  is fit to the treshold point (cross section=0) as well as these two points
*  in attempt to smoothly tail off toward zero.  In some rare cases, however,
*  the cross section is so small and the slope so steep that the parabola can
*  become negative.  In these cases, the differential cross section is set=0.
*
*
*  NEW FOR NSGK:
*  These cross sections do not need to be normalised by the total cross section
*  since these values are the raw double differential cross sections.
*  Interpolation is the same as the above KN setup, for now.
*
*     Revision History:-
*     ================
 
*     2.08     J. Klein     First version
*              Mike Lay     Bug fix index renamed to index1 (conflict with
*                           intrinsic function.
*              N. West      Fix bug reported by Guy Jonkmans re. ICOSINE array.
*     2.09     J. Klein	    Changed logic and added error trap for CC case and
*			    illegal interaction type
*			    Added re-scaling by total cross section.
*              N. West      CCXB,CCXS,NCXB,NXCS banks now all floating point.
*     4.01     S. Brice     Bug fix - NINT added to real pointer increments
*     4.01     I. Lawson    Add capability to read in and use the NSGK cross sections.
*                           Adjust NC threshold from 2.226 to 2.225.
*     4.02     I. Lawson    Fix interpolation for sparse tables.
*     4.02     I. lawson    Added ability to read in the NSAPMGK cross sections.
*                           Fix extrapolation to low energy neutrinos.
 
      IMPLICIT NONE
 
      INCLUDE 'mt.inc'
      INCLUDE 'bank_mctk.inc'
      INCLUDE 'id_particles.inc'
      INCLUDE 'id_en_types.inc'
      INCLUDE 'titles_ncxs.inc'
      INCLUDE 'titles_ncxb.inc'
      INCLUDE 'titles_ccxs.inc'
      INCLUDE 'titles_ccxb.inc'
      INCLUDE 'titles_ccxt.inc'
      INCLUDE 'titles_ncxt.inc'
      INCLUDE 'constants.inc'
      INCLUDE 'zunit.inc'
      INCLUDE 'z.inc'
      INCLUDE 'titles_msw.inc'
 
*     Argument Declarations:-
*     =====================
 
      REAL              ENU,COS_OUT,EOUT
      INTEGER		XSTYPE
 
      REAL ENCCL, ENNCL,ENBCL
*     Local Variable Declarations:-
*     ===========================
 
	REAL E(70),Y(70),COS_ANGLE(12),E_BRACK(2),EOUT_MAX
	REAL DIFXS_E(70),DIFXS_C(12),DIFXS(2),XSEC,SLOPE,XSEC_TOT(2)
	REAL SLOPE1,SLOPE2,THRESH,XS_DIFF
*       Dimension ICOSINE as 13, see note near end of routine.
	INTEGER ICOSINE(13),IBRACK(2),IENERGY,I_OLD,NUTYPE
	INTEGER JENU,JOUT,JCOS,NENU,JJ,LL,L,INDEX1,J
	INTEGER ICNT,ITOT,ITOTOFF,INEXT_ENERGY,CINDEX
*ENDHEADER
 
 
	NUTYPE = IQ(LMCTKP+KMCTK_IDP)
*
*  First find indices for energies that bracket ENU, and save
* energy values
*
        ENCCL = 2.0
        ENNCL = 2.4
        ENBCL = 4.2
 
* Adjust thresholds for NSGK tables, since they are lower
        IF(C_SCTN.EQ.2.OR.C_SCTN.EQ.5) ENNCL = 2.3
 
	IF (XSTYPE.EQ.IDECC) THEN
	  ITOT=LDCCXT
	  IF ((NUTYPE.EQ.IDP_NU_E).AND.(ENU.GE.ENCCL)) THEN
	   IENERGY=LDCCXS+1
	   ITOTOFF=3
	   THRESH=2.225+AMASSP-AMASSN+AMASSE
	   DO WHILE (RCONS(IENERGY).LE.ENU)
	     I_OLD=IENERGY
	     IENERGY=IENERGY+NINT(RCONS((IENERGY-1)+3)) !pointer to next energy value
c             write (iqprnt,10) Enu,
c     +          rcons(i_old),rcons(ienergy)
c10           format (f6.1, 2x, f6.1, 2x, f6.1)

	   ENDDO
	  ELSEIF ((NUTYPE.EQ.IDP_NU_E_BAR).AND.(ENU.GE.ENBCL)) THEN
	   IENERGY=LDCCXB+1
	   ITOTOFF=1
	   THRESH=2.225+AMASSN-AMASSP+AMASSE
	   DO WHILE (RCONS(IENERGY).LE.ENU)
	    I_OLD=IENERGY
	    IENERGY=IENERGY+NINT(RCONS((IENERGY-1)+3)) !pointer to next energy value
	   ENDDO
	  ELSEIF ((NUTYPE.EQ.IDP_NU_E).AND.(ENU.LT.ENCCL)) THEN
	   I_OLD=LDCCXS+1                 !for extrapolation
	   IENERGY=I_OLD+NINT(RCONS((I_OLD-1)+3))
	   ITOTOFF=3
	   THRESH=2.225+AMASSP-AMASSN+AMASSE
	  ELSEIF ((NUTYPE.EQ.IDP_NU_E_BAR).AND.(ENU.LT.ENBCL)) THEN
	   I_OLD=LDCCXB+1                 !for extrapolation
	   IENERGY=I_OLD+NINT(RCONS((I_OLD-1)+3))
	   ITOTOFF=1
	   THRESH=2.225+AMASSN-AMASSP+AMASSE
	  ELSE
	    WRITE(IQPRNT,*)'F_XS_DIFF:WARNING: NEUTRINO TYPE',NUTYPE
            WRITE(IQPRNT,*)'INCONSISTENT WITH CC INTERACTION.'
            WRITE(IQPRNT,*)'DIFFERENTIAL CROSS SECTION WILL BE SET=0.'
	    WRITE(IQLOG,*)'F_XS_DIFF:WARNING: NEUTRINO TYPE',NUTYPE
            WRITE(IQLOG,*)'INCONSISTENT WITH CC INTERACTION.'
            WRITE(IQLOG,*)'DIFFERENTIAL CROSS SECTION WILL BE SET=0.'
	    F_XS_DIFF = 0.
	    GO TO 50
	  ENDIF
	ELSEIF (XSTYPE.EQ.IDENC) THEN !neutral current cross sections needed
	 ITOT=LDNCXT
	 IF (ENU.GE.ENNCL) THEN
	  THRESH=2.225
	  IF ((NUTYPE.EQ.IDP_NU_E).OR.(NUTYPE.EQ.IDP_NU_MU).OR.
     +		(NUTYPE.EQ.IDP_NU_TAU)) THEN
	   IENERGY=LDNCXS+1
	   ITOTOFF=3
	  ELSE
	   IENERGY=LDNCXB+1
	   ITOTOFF=1
	  ENDIF
	  DO WHILE (RCONS(IENERGY).LE.ENU)
	   I_OLD=IENERGY
	   IENERGY=IENERGY+NINT(RCONS((IENERGY-1)+3)) !pointer to next energy value
	  ENDDO
	 ELSE
	  THRESH=2.225
	  IF ((NUTYPE.EQ.IDP_NU_E).OR.(NUTYPE.EQ.IDP_NU_MU).OR.
     +		(NUTYPE.EQ.IDP_NU_TAU)) THEN
	   I_OLD=LDNCXS+1                 !for extrapolation
	   ITOTOFF=3
	  ELSE
	   I_OLD=LDNCXB+1                 !for extrapolation
	   ITOTOFF=1
	  ENDIF
	  IENERGY=I_OLD+NINT(RCONS((I_OLD-1)+3))
	 ENDIF
	ELSE !Illegal cross section type
	    WRITE(IQPRNT,*)'F_XS_DIFF:WARNING: INTERACTION TYPE',XSTYPE
            WRITE(IQPRNT,*)'INCONSISTENT WITH CC OR NC INTERACTION.'
            WRITE(IQPRNT,*)'DIFFERENTIAL CROSS SECTION WILL BE SET=0.'
	    WRITE(IQLOG,*)'F_XS_DIFF:WARNING: INTERACTION TYPE',XSTYPE
            WRITE(IQLOG,*)'INCONSISTENT WITH CC OR NC INTERACTION.'
            WRITE(IQLOG,*)'DIFFERENTIAL CROSS SECTION WILL BE SET=0.'
	    F_XS_DIFF = 0.
	    GO TO 50
	ENDIF
	IBRACK(1)=I_OLD     	     !index of low bracket
	E_BRACK(1)=RCONS(IBRACK(1))  !and energy
	IBRACK(2)=IENERGY            !index of high bracket
	E_BRACK(2)=RCONS(IBRACK(2))  !and energy
*
*  Loop over bracketing energies
	DO JENU=1,2
	 ICOSINE(1)=3    !Index of cosine information offset from energy
 
*  Loop over  cosine values
         CINDEX = 12
	 DO JCOS=1,CINDEX
	  COS_ANGLE(JCOS)=RCONS(IBRACK(JENU)+ICOSINE(JCOS))
	  NENU=RCONS(IBRACK(JENU)+ICOSINE(JCOS)+1)  !number of energy values
	  EOUT_MAX=RCONS(IBRACK(JENU)+ICOSINE(JCOS)+2)
	  IF (EOUT.GT.EOUT_MAX) THEN
	    DIFXS_C(JCOS)=0.0
	  ELSE
*
*  Now store energies and cross section in arrays using Kubodera's
* format.
*
  	   JJ=NENU/8
	   LL=MOD(NENU,8)
	   INDEX1=IBRACK(JENU)+ICOSINE(JCOS)+2 !holds place in RCONS
	   ICNT=1       !holds place in energy and cross section arrays
	   DO J=1,JJ
	    DO JOUT=1,8
	     E(ICNT)=RCONS(INDEX1+JOUT)
	     DIFXS_E(ICNT)=RCONS(INDEX1+8+JOUT)
	     ICNT=ICNT+1
	    ENDDO
 
	    INDEX1=INDEX1+16
           ENDDO
           IF (LL.GT.0) THEN
	    DO L=1,LL
	     E(ICNT)=RCONS(INDEX1+L)
	     DIFXS_E(ICNT)=RCONS(INDEX1+LL+L)
	     ICNT=ICNT+1
	    ENDDO
	   ENDIF
*
*   Spline interpolate to find cross section at given EOUT
*  and current value of COS_ANGLE and bracketing energy E_BRACK(JENU)
          CALL SPLINE(E,DIFXS_E,(ICNT-1),1.E30,1.E30,Y)
          CALL SPLINT(E,DIFXS_E,Y,(ICNT-1),EOUT,XSEC)
          IF(XSEC.LE.0.0) XSEC = 1.0E-15   !set to small number and continue on
	  DIFXS_C(JCOS)=XSEC               !cross section at each cos value
	 ENDIF
*         Note: Last time round loop the next line updates ICOSINE(13).
*         Rather than wasting time protecting this with an IF, just extend
*         array to 13, even though ICOSINE(13) is never used.
	 ICOSINE(JCOS+1)=ICOSINE(JCOS)+3+
     +		2*RCONS(IBRACK(JENU)+ICOSINE(JCOS)+1) !increment pointer
	ENDDO
*
*   Now spline interpolate in cosine
	 CALL SPLINE(COS_ANGLE,DIFXS_C,CINDEX,1.E30,1.E30,Y)
	 CALL SPLINT(COS_ANGLE,DIFXS_C,Y,CINDEX,COS_OUT,XSEC)
	 DIFXS(JENU)=XSEC
	ENDDO
*   NOTE FOR KN CROSS SECTIONS
*   Finally, scale differential values by total cross sections and linearly
*   interplote between bracketing energies to get  value of differential cross
*   section at input values of ENU,COS_OUT, and EOUT
*
*   NOTE FOR NSGK CROSS SECTIONS
*   These cross sections do not need to be scaled by the total cross sections,
*   since they are the pure double differential cross sections. But one must
*   still interpolate to get the differential cross sections at the input values.
 
	 INEXT_ENERGY=1
         XSEC_TOT(1) = 1.0
         XSEC_TOT(2) = 1.0
 
         IF(.NOT.(C_SCTN.EQ.2.OR.C_SCTN.EQ.5)) THEN
	 DO WHILE(RCONS(ITOT+5*INEXT_ENERGY-4).LE.E_BRACK(1))
	   XSEC_TOT(1)=RCONS(ITOT+5*INEXT_ENERGY-ITOTOFF)
	   INEXT_ENERGY=INEXT_ENERGY+1
	 ENDDO
	 XSEC_TOT(2)=RCONS(ITOT+5*INEXT_ENERGY-ITOTOFF)
	 DIFXS(1)=DIFXS(1)*XSEC_TOT(1)
	 DIFXS(2)=DIFXS(2)*XSEC_TOT(2)
         ENDIF
 
	 IF (ENU.GE.E_BRACK(1))  THEN
	   SLOPE = (DIFXS(1)-DIFXS(2))/(E_BRACK(1)-E_BRACK(2))
	   XS_DIFF=(DIFXS(2)+(ENU-E_BRACK(2))*SLOPE)
	   IF (XS_DIFF.GE.0.) THEN
	     F_XS_DIFF=XS_DIFF
	   ELSE                 !we are already zero
	     F_XS_DIFF=0.
	   ENDIF
	 ELSE            !Forced to extrapolate to low energy
	   SLOPE1 = ((ENU-THRESH)*(ENU-E_BRACK(2)))/
     +              ((E_BRACK(1)-THRESH)*(E_BRACK(1)-E_BRACK(2)))
	   SLOPE2 = ((ENU-THRESH)*(ENU-E_BRACK(1)))/
     +              ((E_BRACK(2)-THRESH)*(E_BRACK(2)-E_BRACK(1)))
	   XS_DIFF=DIFXS(1)*SLOPE1 + DIFXS(2)*SLOPE2
 
	   IF (XS_DIFF.GE.0.) THEN
	     F_XS_DIFF=XS_DIFF
	   ELSE                 !we are already zero
	     F_XS_DIFF=0.
	   ENDIF
	 ENDIF
 
50      CONTINUE
	RETURN
	END
*ENDFILE MEMBER=F_XS_DIFF
