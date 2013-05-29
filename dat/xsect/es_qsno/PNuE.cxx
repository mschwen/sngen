////////////////////////////////////////////////////////////////////////
// PNuE
//
// Neutrino - Electron Scattering Cross Section Class
//
// Neutrino Energy, E in MeV
// Electron Energy, T in MeV
//
// Cross sections in 10^42 cm^2
// 
// 17-JAN-2003 modified by Yasuo Takeuchi implement radiative correction
//      o J.N. Bahcall et al PRD 51, 6146-6158 (1995); astro-ph/9502003
//      o modify defaults(), dSigmadT(), Sigma(), Read(), ~Pnue() 
//      o modify SetReaction() to print out mode when it is changed.
//      o clean up PrintReaction() 
// 25-FEB-2003 modified by Yasuo Takeuchi 
//      o kappa and pnc are fixed when fRadiativeCorrection = 2
//      o parameters adjusted with snoman code/get_esdif_rad.for
//      o default = table always
// 03-APR-2003 modified by Yasuo Takeuchi 
//      o replace "double_t" with "Double_t"
// 26-JUN-2003 modified by Yasuo Takeuchi 
//      o fix bug in "fpz" in dSigmadT() 
////////////////////////////////////////////////////////////////////////

#include "PNuE.h"
#include "QNote.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
  
ClassImp(PNuE)
  
  
PNuE::PNuE(Int_t initialize)
{
  // Default Constructor
  
  Defaults();
  fInitialize = initialize;
  if(fInitialize == 1){
    Init();
  }
}
  
PNuE::PNuE(Char_t *rstr)
{
  // Construct a PNuE for the reaction specified in rstr

  Defaults();
  SetReaction(rstr);
  Init();

}

PNuE::~PNuE()
{
    // destructor

    // add by y.t. 16-JAN-2003
    if (fTableTot_e  != NULL) delete fTableTot_e;
    if (fTableTot_mu != NULL) delete fTableTot_mu;
    if (fTableDif_e  != NULL) delete fTableDif_e;
    if (fTableDif_mu != NULL) delete fTableDif_mu;

}
  

void PNuE::Defaults()
{
  // load default parameters
  // M. Chen changed fhbarc, fhbarc2 to RPP 2000, pg. 73
  
  BaseDefaults();

  fDebug = 0;
  //fReaction = nue; // default is nu_e + e -> nu_e + e
  SetReaction("nue");
  fTmin     = 5.0; // minimum allowed recoil energy
 
  fXsSize = 0;
  fEmin   = 0.1;
  fEmax   = 50.0;
  fhbarc  = 197.32696;  // MeV-fm

  //------------------------------------------------------------------------ 
  // these parameters should be adjusted with SNOMAN code/get_esdif_rad.for
  //------------------------------------------------------------------------ 
  me      = 0.510998918;    // MeV/c2; updated to PDG2004 by JH
  fGf     = 1.16637e-5;     // GeV-2; updated to PDG2004 by JH
  fhbarc2 = 0.389379323;    // GeV^2 mb; updated to PDG2004 by JH
  fpi     = 3.1415926535897932; // more digits, added by JH
  falpha  = 0.007297352568; // fine structure; updated to PDG2004 by JH
  fsinthetaW2 = 0.23120;  // Weinburg angle; updated to PDG2004 by JH
  //------------------------------------------------------------------------ 
  
  sigmaoverme = 2*fGf*fGf*me*1.0e9*fhbarc2/fpi; // 10^-42 cm^2/MeV

  // add by y.t. 16-JAN-2003  flag to use radiative correction or not
  fRadiativeCorrection = 4; //  ==1     previous routine
                            //  ==2     new routine w/o  rad. cor. 
                            //  ==3     new routine with rad. cor. w/o table
                            //  ==4     new routine with rad. cor. with table
                            //  others  new routine with rad. cor. with or  
                            //          w/o table (fTableTeMin/fTableTeMax)

  fTableTeMin = 0.05;       // don't use table from 0 upto 0.05*Temax
  fTableTeMax = 0.90;       // don't use table from 0.90*Temax upto Temax

  fTableTot_e  = NULL;
  fTableTot_mu = NULL;
  fTableDif_e  = NULL;
  fTableDif_mu = NULL;

  /*
  fColor = kBlue;
  fLineWidth = 1;
  fLineStyle = 1;
  fColor = kRed;
  fMarkerStyle = 21;
  fColor   = 19;
  */
}
  
void PNuE::Init()
{
  // initialize NuE object

  fNote = new QNote(fNoteDepth);
  if(fVerbosity >= 1){
    Char_t message[132];
    sprintf(message,"creating PNuE object, reaction = %s",fReactionString);
    fNote->Note(message);
  }

  BaseInit();

  //PrintReaction();
  CalcG();

  // add by y.t. 16-JAN-2003  
  // read cross section tables with radiative correction
  Read();

}

void PNuE::SetReaction(Char_t *rstr)
{
  // Set the reaction type
  // valid types are:  nue     --  nu_e + e       -> nu_e + e
  //                   nuebar  --  anti nu_e + e  -> anti nu_e + e
  //                   numu    --  nu_mu + e      -> nu_mu + e
  //                   numubar --  anti nu_mu + e -> anti nu_mu + e

  if( !strcmp(rstr,"nue")){
    fReaction = nue;
    sprintf(fReactionString,"nu_e+e-->nu_e+e");
  } else if( !strcmp(rstr,"nuebar")){
    fReaction = nuebar;
    sprintf(fReactionString,"!nu_e+e-->!nu_e+e");
  } else if( !strcmp(rstr,"numu")){
    fReaction = numu;
    sprintf(fReactionString,"nu_mu+e-->nu_mu+e");
  } else if( !strcmp(rstr,"numubar")){
    fReaction = numubar;
    sprintf(fReactionString,"!nu_mu+e-->!nu_mu+e");
  } else {
    printf("  unknown reaction type: %s\n",rstr);
    fReaction = nue;
    sprintf(fReactionString,"nu_e+e-->nu_e+e");
  }

  fNote = new QNote(fNoteDepth);
  if(fVerbosity >= 1){
    Char_t message[132];
    sprintf(message,"setting reaction to %s", fReactionString);
    fNote->Note(message);
  }

  CalcG();

}
  
void PNuE::PrintReaction()
{
  // Set the reaction type
  // valid types are:  nue     --  nu_e + e       -> nu_e + e
  //                   nuebar  --  anti nu_e + e  -> anti nu_e + e
  //                   numu    --  nu_mu + e      -> nu_mu + e
  //                   numubar --  anti nu_mu + e -> anti nu_mu + e

  // fReactionString is set in SetReaction()  by y.t. 14-JAN-2003
  //if(        fReaction == nue){
  //  sprintf(fReactionString,"nu_e+e-->nu_e+e");
  //  //printf("   nu_e + e    -->  nu_e + e\n");
  //} else if( fReaction == nuebar){
  //  sprintf(fReactionString,"!nu_e+e-->!nu_e+e");
  //  //printf("   !nu_e + e   -->  !nu_e + e\n");
  //} else if( fReaction == numu){
  //  sprintf(fReactionString,"nu_mu+e-->nu_mu+e");
  //  //printf("   nu_mu + e   -->  nu_mu + e\n");
  //} else if( fReaction == numubar){
  //  sprintf(fReactionString,"!nu_mu+e-->!nu_mu+e");
  //  //printf("   !nu_mu + e  -->  !nu_mu + e\n");
  //}

  printf("   %s", fReactionString);

}


void PNuE::CalcG()
{
  // calculate the gL and gR for 
  // the reaction type

  if ( fReaction == nue)           {
  
    fgL =  0.5 + fsinthetaW2;
    fgR =  fsinthetaW2;
  
  } else if (fReaction == nuebar)  {

    fgL =   fsinthetaW2;
    fgR =   0.5 + fsinthetaW2;

  } else if (fReaction == numu)    {

    fgL =   -0.5 + fsinthetaW2;
    fgR =    fsinthetaW2;

  } else if (fReaction == numubar) {
  
    fgL =   fsinthetaW2;
    fgR =   -0.5 + fsinthetaW2;

  }
}

//--------------------------------------------------------------
// Read() part rewrite by y.t. 16-JAN-2003  

void PNuE::Read()
{
    // Read from default file
    Read(SIGMA_TOT_FILE, SIGMA_DIF_FILE);
}

void PNuE::Read(Char_t *filename1, Char_t *filename2)
{
    // Read the cross section files with radiative correction
    // Use the QDPath default path to find file.

    FILE *fp1, *fp2;
  
    if(fVerbosity >=  1){
	printf("PNuE: read table from total:%s  diff:%s \n", 
	       filename1, filename2);
    }
    fp1 = gQDPath->Open(filename1,"r");
    if(fp1 == 0){
	printf("PNuE: can't open total:%s\n", filename1);
	exit(0);
    }
    fp2 = gQDPath->Open(filename2,"r");
    if(fp2 == 0){
	printf("PNuE: can't open diff:%s\n", filename2);
	exit(0);
    }
    Read(fp1, fp2); 
    fclose(fp1);
    fclose(fp2);
}

void PNuE::Read(FILE *fp1, FILE *fp2)
{
    Char_t line[MAX_LINE_BUF];
    Double_t Enu, Te, sigma_e, sigma_mu;
    Int_t k;

    // Initialize table 
    if(fVerbosity >=  1){
	printf("PNuE: reset tables.\n");
    }

    fNDataTot = 0;           
    fEnuStepTot = 0.0;
    fNDataDif = 0;           
    fEnuStepDif = 0.0;
    if (fTableTot_e  != NULL) delete fTableTot_e;
    if (fTableTot_mu != NULL) delete fTableTot_mu;
    if (fTableDif_e  != NULL) delete fTableDif_e;
    if (fTableDif_mu != NULL) delete fTableDif_mu;

    // Read total cross section from fp1
    while (fgets(line, MAX_LINE_BUF, fp1)!= NULL) {
	if (line[0] != '#') {
	    if (fEnuStepTot == 0.0) {
		sscanf(line, "EnuStepTot = %lf\n", &fEnuStepTot);
	    } else if (fNDataTot == 0) {
		sscanf(line, "NDataTot   = %d\n", &fNDataTot);
		fTableTot_e  = new double[fNDataTot];
		fTableTot_mu = new double[fNDataTot];
	    } else {
		sscanf(line, "  %lf  %le  %le\n", &Enu, &sigma_e, &sigma_mu);
		k = ((int)(Enu/fEnuStepTot+0.5) - 1);
		fTableTot_e[k] = sigma_e;
		fTableTot_mu[k] = sigma_mu;

		if (fVerbosity >= 2) {
		    line[strlen(line)-1]= '\0';
		    printf("%d: line=<%s> Enu=%f sigma_e=%e sigma_mu=%e\n",
			   k, line, Enu, sigma_e, sigma_mu);
		}
	    }
	}
    }
	  
    if (fVerbosity >= 1) {
	printf("PNuE: total: fNDataTot=%d  fEnuStepTot=%f  last k=%d\n", 
	       fNDataTot, fEnuStepTot, k);
    }

    // Read differential cross section from fp2
    while (fgets(line, MAX_LINE_BUF, fp2)!= NULL) {
	if (line[0] != '#') {
	    if (fEnuStepDif == 0.0) {
		sscanf(line, "EnuStepDif = %lf\n", &fEnuStepDif);
	    } else if (fNDataDif == 0) {
		sscanf(line, "NDataDif   = %d\n", &fNDataDif);
		fTableDif_e  = new double[fNDataDif*fNDataDif];
		fTableDif_mu = new double[fNDataDif*fNDataDif];
	    } else {
		sscanf(line, "  %lf  %lf  %le  %le\n",
		       &Enu, &Te, &sigma_e, &sigma_mu);
		k = ((int)(Enu/fEnuStepDif+0.5) - 1) * fNDataDif + 
		    ((int)(Te/fEnuStepDif+0.5) - 1);
		fTableDif_e[k]  = sigma_e;
		fTableDif_mu[k] = sigma_mu;

		if (fVerbosity >= 2) {
		    line[strlen(line)-1]= '\0';
		    printf("%d: line=<%s> Enu=%f Te=%f sig_e=%e sig_mu=%e\n",
			   k, line, Enu, Te, sigma_e, sigma_mu);
		}
	    }
	}
    }
	  
    if (fVerbosity >= 1) {
	printf("PNuE: diff: fNDataDif=%d  fEnuStepDif=%f  last k=%d\n", 
	       fNDataDif, fEnuStepDif, k);
    }

}

//-------------------------------------------------------------------

// modified by y.t. 16-JAN-2003
Double_t PNuE::Sigma(Double_t Enu)
{
    // return total cross section in units cm^-42
    // for laboratory neutrino energy Enu

    Double_t   sigma = 0.0;
    Double_t   Temax = Enu - 1.0/(2.0/me + 1.0/Enu);

    if (fRadiativeCorrection == 1) { 

	/////////////////////////////////////////////////////////////
	// previous routine
	/////////////////////////////////////////////////////////////
	//
	// modified by M. Chen to integrate up to the kinematic maximum
	
	Double_t   gL2 = fgL*fgL;
	Double_t   gR2 = fgR*fgR;
	Double_t   gLR = fgL*fgR;

	sigma =    (gL2 + gR2)*Temax
	    - (gR2/Enu + gLR*me/(2.0*Enu*Enu))*Temax*Temax
	    + (gR2/(3.0*Enu*Enu))*Temax*Temax*Temax;

	sigma = sigma*sigmaoverme;

	/* ******************** *\
	   Double_t g = fgL*fgL + (fgR*fgR)/3.0;  
	   sigma = 1.0e42 * (2.0 * fGf*fGf * me * elab * g 
	               * fhbarc2 * 1.0e-6 * 1.0e-27) / fpi;
	\* ******************** */

	if(fDebug){
	    printf("  PNuE::Sigma(%f)\n",Enu); 
	    printf("     gL = %f  gR = %f\n",fgL,fgR);
	    printf("     sigma=%f\n",sigma);
	}

	return sigma;

    } else {

	/////////////////////////////////////////////////////////////
	// radiative correction
	/////////////////////////////////////////////////////////////
	Int_t k = ((int)(Enu/fEnuStepTot) - 1);

	if (fRadiativeCorrection != 3 && fRadiativeCorrection != 2 &&
	    k >= 0 && k < fNDataTot - 1) {

	    // use table
	    Double_t e1 = fEnuStepTot * (double) (k+1);
	    Double_t e2 = fEnuStepTot * (double) (k+2);

	    if ( fReaction == nue) {
		// interpolate by a straight line
		sigma = (fTableTot_e[k+1] - fTableTot_e[k]) / 
		    (e2 - e1) * (Enu - e1) + fTableTot_e[k];
	    } else if (fReaction == numu)    {
		// interpolate by a straight line
		sigma = (fTableTot_mu[k+1] - fTableTot_mu[k]) / 
		    (e2 - e1) * (Enu - e1) + fTableTot_mu[k];
	    } else {
		// not supported yet
		cout << "PNuE: error: " << fReaction; 
		cout << " is not supported!!\n";
		exit(-1);
	    }

	    if (fVerbosity >= 2) {
		printf("e1=%f e2=%f se1=%e se2=%e smu1=%e smu2=%e: Enu=%f sig=%e\n", 
		       e1, e2, fTableTot_e[k], fTableTot_e[k+1],
		       fTableTot_mu[k], fTableTot_mu[k+1],
		       Enu, sigma);
	    }

	    return sigma;

	} else {

	    // calculate now (about 1.3 sec / calculation...)

	    // To avoid log(-1.0e-15) in dSigmadT()
	    Temax = Temax *(1.0 - 1.0e-13);
	    
	    Int_t     istep = 1000;     // should not be hard corded??
	    Double_t  dstep = Temax / (double) istep;
	    Double_t  Te;

	    for (Int_t i = 1; i <= istep; i++) {
		Te = dstep * (double) i;
		sigma = sigma + dSigmadT(Enu, Te) * dstep;
	    } 
	    if(fDebug){
		printf("PNuE::Sigma(%g): Temax=%g sig=%g\n",
		       Enu, Temax, sigma);
	    }
	    return sigma;
	}

    }

    // should not come here
    return -1.0;

}

Double_t PNuE::SigmaLab(Double_t elab)
{
  return Sigma(elab);
}

//--------------------------------------


// modified by y.t. 14-JAN-2003
Double_t PNuE::dSigmadT(Double_t Enu,Double_t Te)
{
    // Differential cross section
    //  d Sigma
    //  -------
    //  d Te
    //
    //  Enu = laboratory neutrino energy        -- units: MeV
    //  Te  = laboratory recoil electron energy -- units: MeV
    //  dSigmadT units:  10^-42 cm^2/MeV
    
    if (fRadiativeCorrection == 1) { 
	
	/////////////////////////////////////////////////////////////
	// previous routine
	/////////////////////////////////////////////////////////////

	// Taken from Bahcall "Neutrino Astrophysics" eqn 8.29 and 8.31
	//
	// Fraser's coding hacked from Bahcall, Kamionkowski and Sirlin
	// is correct, except that the kinematic limit (Compton edge)
	// needs to be added: M. Chen, 14 March 2001
	//
	// radiative corrections are NOT included in this form
	//
	// M. Chen added:
	//   - computation of 17.216e-3 from fundamental constants
	//   - sigmaoverme is calculated in the default constructor
	//   - kinematic limit for Te

        /* ******************** *\
	   Double_t Sigma = 88.083e-46 * 1.0e42; // sigma in units 10^-42
	   
	   Double_t     y = Te/Enu;
	   Double_t   gL2 = fgL*fgL;
	   Double_t   gR2 = fgR*fgR;
	   Double_t   gLR = fgL*fgR;
	   Double_t oneMy = 1.0 - y;
	   

	   Double_t dsigdT = Sigma * ( gL2 + (gR2*oneMy*oneMy) - (gLR*y) );
  
	   return dsigdT;  
	\* ******************** */
 
	// hacked together from eqn A1 of Bahcall, Kamionkowski and Sirlin
	// astro-ph/9502003

	Double_t Temax = Enu - 1.0/(2.0/me + 1.0/Enu);
	
	if (Te > Temax)
	    return 0.0;
    
	Double_t z = Te/Enu;
	Double_t mq = me/Enu;
	
	Double_t   gL2 = fgL*fgL;
	Double_t   gR2 = fgR*fgR;
	Double_t   gLR = fgL*fgR;
	Double_t   oneMz = 1.0 - z;
  
	Double_t dsigdT = sigmaoverme * (gL2 + (gR2*oneMz*oneMz) - gLR*mq*z);
  
	// Double_t dsigdT = 17.216e-3 * ( gL2 + (gR2*oneMz*oneMz) - gLR*mq*z);

	if (fDebug) {
	    printf("Enu=%f(MeV) Te=%f(MeV) dsigma/dT=%g x10-42 cm2/MeV\n",
		   Enu, Te, dsigdT);
	}

	return dsigdT;  
    }
    else {

	/////////////////////////////////////////////////////////////
	// radiative correction
	/////////////////////////////////////////////////////////////

	// check kinematical limit
	Double_t Temax = Enu - 1.0/(2.0/me + 1.0/Enu);

	if (Te < 0.0 || Te > Temax) {
	    return 0.0;
	}

	// check Te == 0.0 (Is this correct?)
	if (Te == 0.0) {
	    if(fVerbosity >= 1){
		printf("PNuE: warning: Te==0.0 ???, use 1eV\n");
	    }
	    Te = 0.000001;    // 1eV
	}

	// check table 
	Int_t k1 = ((int)(Enu/fEnuStepDif) - 1);     
	Int_t k2 = ((int)(Te/fEnuStepDif) - 1);      

	Double_t e1 = fEnuStepDif * (double) (k1+1); 
	Double_t e2 = fEnuStepDif * (double) (k1+2); 
	Double_t t1 = fEnuStepDif * (double) (k2+1);
	Double_t t2 = fEnuStepDif * (double) (k2+2);

	Double_t sig_e1_t1 = 0.0;
	Double_t sig_e1_t2 = 0.0;
	Double_t sig_e2_t1 = 0.0;
	Double_t sig_e2_t2 = 0.0;

	Double_t Te1 = Temax * fTableTeMin;
	Double_t Te2 = Temax * fTableTeMax;

	if (k1 >= 0 && k1 < fNDataDif - 1 &&
	    k2 >= 0 && k2 < fNDataDif - 1 ) {

	    if ( fReaction == nue) {
		sig_e1_t1 = fTableDif_e[k1 * fNDataDif + k2];  
		sig_e1_t2 = fTableDif_e[k1 * fNDataDif + k2+1];  
		sig_e2_t1 = fTableDif_e[(k1+1) * fNDataDif + k2];  
		sig_e2_t2 = fTableDif_e[(k1+1) * fNDataDif + k2+1];  
	    } else if (fReaction == numu)    {
		sig_e1_t1 = fTableDif_mu[k1 * fNDataDif + k2];  
		sig_e1_t2 = fTableDif_mu[k1 * fNDataDif + k2+1];  
		sig_e2_t1 = fTableDif_mu[(k1+1) * fNDataDif + k2];  
		sig_e2_t2 = fTableDif_mu[(k1+1) * fNDataDif + k2+1];  
	    } else {
		// not supported yet
		cout << "PNuE: error: " << fReaction; 
		cout << " is not supported!!\n";
		exit(-1);
	    }

	}

	if (fRadiativeCorrection != 3 && fRadiativeCorrection != 2 &&
	    k1 >= 0 && k1 < fNDataDif - 1 &&
	    k2 >= 0 && k2 < fNDataDif - 1 &&
	    sig_e1_t1 > 0.0 && sig_e1_t2 > 0.0 &&
	    sig_e2_t1 > 0.0 && sig_e2_t2 > 0.0 &&
	    (fRadiativeCorrection == 4 || (Te > Te1 && Te < Te2))) {

	    // use table
	    Double_t r1 = (sig_e1_t2 - sig_e1_t1) / (t2 -t1) 
		* (Te - t1) +  sig_e1_t1;
	    Double_t r2 = (sig_e2_t2 - sig_e2_t1) / (t2 -t1) 
		* (Te - t1) +  sig_e2_t1;
            /** 03-APR-2003 pointed out by Neil, fixed by y.t. 
	    double_t r1 = (sig_e1_t2 - sig_e1_t1) / (t2 -t1) 
		* (Te - t1) +  sig_e1_t1;
	    double_t r2 = (sig_e2_t2 - sig_e2_t1) / (t2 -t1) 
		* (Te - t1) +  sig_e2_t1;
 	    **/
	    return ((r2 - r1) / (e2 - e1) * (Enu - e1) + r1);

	}
	else {
	    // don't use table
	    Double_t E = me + Te;
	    Double_t x = sqrt(1.0 + 2.0 * me / Te);
	    Double_t z = Te / Enu;
	    Double_t el = sqrt(E * E - me * me);
	    Double_t IT = 1.0 / 6.0 * (1.0 / 3.0 + (3.0 - x * x) 
			  * (0.5 * x * log((x + 1.0) / (x - 1.0)) - 1.0));
	 
	    // To avoid log(-1.0e15) in fm,fpz,fpm when Te = Temax 
	    if (1.0-z-me/(E+el) <= 0.0) {
		if(fVerbosity >= 1){
		    printf("PNuE: warning: 1.0-z-me/(E+el) = %e\n",
		       1.0-z-me/(E+el));
		}
		return 0.0;
	    }

	    // just use central values for pnc and kappa (correct ?)
	    Double_t gl, gr, kappa, pnc = 1.0126;

	    if ( fReaction == nue) {
		kappa = 0.9791 + 0.0097 * IT;
		if (fRadiativeCorrection == 2) {
		    pnc = kappa = 1.0;
		}
		gl =  pnc * (0.5 - kappa * fsinthetaW2) - 1.0;
		gr = -pnc * kappa * fsinthetaW2;
	    } else if (fReaction == numu)    {
		kappa = 0.9970 - 0.00037 * IT;
		if (fRadiativeCorrection == 2) {
		    pnc = kappa = 1.0;
		}
		gl =  pnc * (0.5 - kappa * fsinthetaW2);
		gr = -pnc * kappa * fsinthetaW2;

	    } else {
		// not supported yet
		cout << "PNuE: error: " << fReaction; 
		cout << " is not supported!!\n";
		exit(-1);
	    }

	    // turn off radiative correction, if fRadiativeCorrection == 2
	    Double_t fm  = 0.0;
	    Double_t fpz = 0.0;
	    Double_t fpm = 0.0;

	    if (fRadiativeCorrection != 2) {
		fm = (E/el * log((E+el)/me)-1.0) *
		    (2.0*log(1.0-z-me/(E+el))
		     - log(1.0-z) - 0.5*log(z) - 5.0/12.0) + 0.5 *
		    (fL(z) - fL(el/E)) - 0.5 * log(1.0-z)*log(1.0-z) -
		    (11.0/12.0 + 0.5 * z) * log(1.0-z)
		    + z * (log(z) + 0.5 * log(2.0*Enu/me))
		    - (31.0/18.0 + 1.0/12.0 * log(z))*el/E
		    - 11.0/12.0 * z + z*z/24.0;
	    
		fpz = (E/el * log((E+el)/me)-1.0) * 
		    ((1.0-z)*(1.0-z) * (2.0*log(1.0-z-me/(E+el))  
		    - log(1.0-z) - 0.5*log(z) - 2.0/3.0) - 0.5 * 
		     (z*z * log(z) + 1.0 - z))
		    - 0.5 * (1.0-z)*(1.0-z) * (log(1.0-z)*log(1.0-z) + el/E * 
		     (fL(1.0-z) - log(z) * log(1.0-z))) 
		    + log(1.0-z) * (z*z/2.0*log(z) + (1.0-z)/3.0*
		     (2.0*z-1.0/2.0))
		    - 0.5 * z*z * fL(1.0-z) - z*(1.0-2.0*z)/3.0 * 
		    log(z) - z*(1.0-z)/6.0 - el/E/12.0*(log(z) + 
		     (1.0-z)*(115.0-109.0*z)/6.0);
	    
		fpm = (E/el * log((E+el)/me)-1.0) * 2.0 *
		    log(1.0-z-me/(E+el));
	    }

	    Double_t dsigma_dT = 2.0 * fGf*fGf * me / fpi * (
		gl * gl * (1.0 + falpha / fpi * fm)
		+ gr * gr * ((1.0 - z) * (1.0 - z) + falpha / fpi * fpz)
		- gr * gl * me / Enu * z * (1.0 + falpha / fpi * fpm));

	    // 26-JUN-2003 bug fixed by y.t.     
	    // fpz already contains (1.0-z)**2 factor
	    // + gr * gr * (1.0 - z) * (1.0 - z) * (1.0 + falpha / fpi * fpz)

	    if (fDebug) {
		printf("Enu=%f(MeV) Te=%g(MeV) dsigma/dT=%g x10-42 cm2/MeV\n",
		       Enu, Te, dsigma_dT * fhbarc2 * 1.0e9);
		printf("Temax=%g fm=%g  fpz=%g  fpm=%g\n", 
		       Temax, fm, fpz, fpm);
		printf("x=%g IT=%g  gl=%g  gr=%g\n", 
		       x, IT, gl, gr);	
		printf("(E+el)/me=%g z=%g 1.0-z-me/(E+el)=%g\n", 
		       (E+el)/me, z, 1.0-z-me/(E+el));	

	    }
	    return   dsigma_dT * fhbarc2 * 1.0e9; 
	}

    }

    // should not come here
    return -1.0;
}


void PNuE::DrawdSigmadT(Double_t Enu,Char_t *options)
{
  
  // Draw the xs

  Int_t i;
  Float_t x[MAXPOINTS],y[MAXPOINTS];

  Int_t numpoints = 500;
  
  Double_t Emax = Enu;
  Double_t Emin =   0.0;

  Double_t xwid = Emax - Emin;
  Double_t xstep = xwid/(Double_t)numpoints;

  for(i=0;i<=numpoints;i++){
    x[i] = Emin + (Double_t)i * xstep;
    y[i] = dSigmadT(Enu,x[i]);
  }  

  //for(i=0;i<fXsSize;i++){
  //  x[i] = fXs[0][i];
  //  y[i] = fXs[1][i];
  //}

  if(fgr != 0){
    delete fgr;
  }
  fgr = new TGraph(numpoints+1,x,y);
  fgr->SetFillColor(fColor);
  fgr->SetLineColor(fColor);
  fgr->SetLineWidth(fLineWidth);
  fgr->SetMarkerColor(fColor);
  fgr->SetMarkerStyle(fMarkerStyle);
  Char_t title[132];
  sprintf(title,"%s(Enu=%0.3gMeV)\n",fReactionString,Enu);
  fgr->SetTitle(title);
  fgr->Draw(options);
  if( FindOption(options,'a') || FindOption(options,'A')){
    gPad->Update();
    fgr->GetHistogram()->SetXTitle("Te(MeV)");
    fgr->GetHistogram()->SetYTitle("dSigma/dTe(10^-42 cm^2/MeV)");
  }
}

void PNuE::DrawdSigmadT(Double_t Enu)
{
  DrawdSigmadT(Enu,"AL");
}



Double_t PNuE::IntegraldSigmadT(Double_t Enu,Double_t T1,Double_t T2)
{
  // Integrate dSigma/dT of Enu from T1 to T2
  //
  // the code below could be replaced someday with
  // analytic integration, as opposed to numerical
  // M. Chen, 14 March 2001
 
  Int_t i;
  Int_t  nbins = 100;
  Double_t integral;
  Float_t x1 = T1;
  Float_t x2 = T2;
  Float_t T,dsig;
  Float_t Tstep = (x2-x1)/(Float_t)nbins;

  //TF1 *f = new TF1("yield",(Double_t)Eval2,x1,x2,0);
  //Double_t integral = f->Integral(x1,x2);
  //printf("  Integral = %f\n",integral);
  
  TH1F *h = new TH1F("h1","h1",nbins,x1,x2);
  for(i=0, T=x1;i<nbins;i++,T+= Tstep){
    dsig = dSigmadT(Enu,T);
    h->Fill(T + (Tstep/2.0),dsig * Tstep);
  }
  integral = h->Integral();
  //printf("  E1=%f E2=%f Estep = %f\n",E1,E2,Estep);
  printf("  Integral = %f\n",integral);
  //h->Draw();
  delete h;

  return(integral);
}


//--------------------------------------


Double_t PNuE::dSigmadCosTh(Double_t Enu,Double_t CosTh)
{
  // Differential cross section
  //  d Sigma
  //  -------
  //  d CosTh
  //
  //  Enu = laboratory neutrino energy        -- units: MeV
  //  CosTh = laboratory recoil electron direction cosine
  //  dSigmadT units:  10^-42 cm^2/MeV
  //
  // Taken from Bahcall "Neutrino Astrophysics" eqn 8.40
  //
  // modified: M. Chen 14 March 2001
  // correcting John Bahcall's bad Te in units of me formulae
 



  Double_t mu = CosTh;
  Double_t mu2 = mu*mu;
  Double_t Enu2 = Enu*Enu;
  Double_t meEnu = me + Enu;

  // first see if scattering angle is greater than 90 degrees
  // (not allowed)
  if(mu < 0.0){
    return 0.0;
  }



  Double_t Denom = meEnu*meEnu - Enu2*mu2;

  Double_t T  = (2.0 * Enu2 * mu2 * me) / Denom;

  Double_t dsigdT = dSigmadT(Enu,T);
  
  Double_t dSigmadCosTh = dsigdT *me*4.0*meEnu*meEnu*Enu2*mu/(Denom*Denom);

  if(fDebug){
    printf("  dSigmadCosTh(%g,%g)\n",Enu,CosTh);
    printf("    mu=%g mu2=%g Enu2=%g\n",mu,mu2,Enu2);
    printf("    Denom = %g\n",Denom);
    printf("    T = %g dSigmadT = %g\n",T,dsigdT);
    printf("    dSigmadCosTh = %g\n",dSigmadCosTh);
  }
  
  return dSigmadCosTh;
}


void PNuE::DrawdSigmadCosTh(Double_t Enu,Char_t *options)
{
  
  // Draw the xs

  Int_t i;
  Float_t x[MAXPOINTS],y[MAXPOINTS];

  Int_t numpoints = 500;

  Double_t xwid = 2.0;
  Double_t xstep = xwid/(Double_t)numpoints;

  for(i=0;i<=numpoints;i++){
    x[i] = -1.0 + (Double_t)i * xstep;
    y[i] = dSigmadCosTh(Enu,x[i]);
  }  

  //for(i=0;i<fXsSize;i++){
  //  x[i] = fXs[0][i];
  //  y[i] = fXs[1][i];
  //}

  if(fgr != 0){
    delete fgr;
  }
  fgr = new TGraph(numpoints+1,x,y);
  fgr->SetFillColor(fColor);
  fgr->SetLineColor(fColor);
  fgr->SetLineWidth(fLineWidth);
  fgr->SetMarkerColor(fColor);
  fgr->SetMarkerStyle(fMarkerStyle);
  Char_t title[132];
  sprintf(title,"%s(Enu=%0.3gMeV)\n",fReactionString,Enu);
  fgr->SetTitle(title);
  fgr->Draw(options);
  if( FindOption(options,'a') || FindOption(options,'A')){
    gPad->Update();
    fgr->GetHistogram()->SetXTitle("cos(Th)");
    fgr->GetHistogram()->SetYTitle("dSigma/dCosTh(10^-42 cm^2)");
  }
}

void PNuE::DrawdSigmadCosTh(Double_t Enu)
{
  DrawdSigmadCosTh(Enu,"AL");
}


//----------------------------------------------

Double_t PNuE::dSigmadOmega(Double_t e,Double_t theta,Double_t phi)
{
  // return the angular differential cross section in units
  // cm^2/str
  
  return 0;
}


//---------------- private routine ----------------------------

// fL(): add by y.t. 14-JAN-2003
Double_t PNuE::fL(Double_t x)
{
    int istep = 1000;    // should not be hard corded??

    double dstep = x / (double) istep;
    double t;
    double sum = 0.0;

    for (int i = 1; i <= istep; i++) {
         t = dstep * (double) i;
         sum = sum + log(fabs(1.0 - t))/t * dstep;
    } 

    return sum;
}


/**********************************************************************:
 **** comment out by y.t. 16-JAN-2003 these routine seems not using...
 **** so, rewrite to read radiative correction tables

void PNuE::Read(char *filename)
{
  // Read xs from a file
  // File format:
  //  #
  //  # comments
  //  #
  //  # 
  //  x1  y1
  //  x2  y2
  //  x3  y3
  //   ...
  //
  //  The xs is assumed to be monotomically increasing
  //  in x.
  

  FILE *fp;
  Char_t line[256];
  
  
  fp = fopen(filename,"r");
  if(fp == 0){
    printf("  can't open %s\n",filename);
    exit(0);
  }
  
  while (GetLine(fp,line,'#')){
    //if(fDebug){
    //  printf("line=<%s>\n",line);
    //}
    if(line[0] != '\0'){
      sscanf(line,"%lf %lf",
         &(fXs[0][fXsSize]),&(fXs[1][fXsSize]));
      fXsSize ++;
    }
    if(fXsSize >= MAXPOINTS){
      printf("  PNuE:: Maximum Number of xs points reached.\n");
      break;
    }
  }
  fEmin = fXs[0][0];
  fEmax = fXs[0][fXsSize-1];

  printf("  %d point xs read from %s\n",fXsSize,fFileName);
  printf("  range [%f,%f]\n",fEmin,fEmax);
  
  fclose(fp);
  
  //fSpec = new TF1(filename,Get2,fEmin,fEmax,0);

}
**************************************************************************/
  

/*
Double_t PNuE::Eval(Double_t x)
{
  // Return the xs value at x.
  // If x is outside the range of the xs, return
  // the closest xs value.
  // Otherwise, do a linear interpolation between points.
  
  Int_t   i;
  Int_t   ilow  = 0;
  Int_t   ihigh = 0;
  Double_t xlow, xhigh;
  Double_t ylow, yhigh;
  Double_t xwid, ywid;
  Double_t xdif, ydif;
  Double_t y = 0;

  if( x <= fEmin){
    return fXs[1][0];
  } else if ( x >= fEmax) {
    return fXs[1][fXsSize-1];
  } else {
    for(i=0; i < fXsSize; i++){
      if(fXs[0][i] > x){
        ihigh = i;
        ilow  = i - 1;
        break;
      }
    }
    xlow  = fXs[0][ilow];
    xhigh = fXs[0][ihigh];
    ylow  = fXs[1][ilow];
    yhigh = fXs[1][ihigh];
  
    xwid  = xhigh - xlow;
    ywid  = yhigh - ylow;
  
    xdif  = x - xlow;
    ydif  = (xdif/xwid) * ywid;
    y     = ylow + ydif;  

    if(fDebug){
      printf("  PNuE::Get(%f)\n",x);
      printf("    xlow  = %f ylow  = %f\n",xlow,ylow);
      printf("    xhigh = %f yhigh = %f\n",xhigh,yhigh);
      printf("    x     = %f y     = %f\n",x,y);
    }
  }
  return y;
}
*/  

/*
void PNuE::Draw(Char_t *options)
{
  
  // Draw the total cross section (units 10^-42 cm^2)

  Int_t i;
  Float_t x[MAXPOINTS],y[MAXPOINTS];

  Int_t numpoints = 500;
  
  Double_t xwid = fEmax - fEmin;
  Double_t xstep = xwid/(Double_t)numpoints;

  for(i=0;i<=numpoints;i++){
    x[i] = fEmin + (Double_t)i * xstep;
    y[i] = Sigma(x[i]);
  }  

  //for(i=0;i<fXsSize;i++){
  //  x[i] = fXs[0][i];
  //  y[i] = fXs[1][i];
  //}

  if(fgr != 0){
    delete fgr;
  }
  fgr = new TGraph(numpoints+1,x,y);
  fgr->SetFillColor(fColor);
  fgr->SetLineColor(fColor);
  fgr->SetLineStyle(fLineStyle);
  fgr->SetLineWidth(fLineWidth);
  fgr->SetMarkerColor(fColor);
  fgr->SetMarkerStyle(fMarkerStyle);
  fgr->SetTitle(fReactionString);
  fgr->Draw(options);
  if( FindOption(options,'a') || FindOption(options,'A')){
    gPad->Update();
    fgr->GetHistogram()->SetXTitle("E_{#nu} (MeV)");
    fgr->GetHistogram()->SetYTitle("#sigma (cm^{2}) #times 10^{42}");
  }
}

void PNuE::Draw()
{
  Draw("AL");
}
*/


//-------------------------------------------------------------
// private routines

  
/*

Int_t PNuE::GetLine(FILE *fp,char *line,char comchar)
{
  // Get an input line from *fp.  Only return everything up to the
  // comchar.
    char c;
    int i;
  
    for(i=0;((c=getc(fp))!='\n')&&(c!=EOF)&&(c!=comchar);i++){
        line[i]=c;
    }
    if(c==comchar)
        while(((c=getc(fp))!='\n')&&(c!=EOF))
            ;
    line[i]='\0';
    if(c==EOF)
        return(0);
    else 
        return(1);
}
   
  

Int_t PNuE::FindOption(Char_t *options,Char_t option)
{
  // Look for a particular plotting option
  // return 1 if found.
  
  Char_t c;
  Int_t i;
  for(i=0; (c=options[i]) !='\0';i++){
    if(fDebug){
      printf("  FindOptions(<%s>,%c) c=%c\n",options,option,c);
    }
    if( c == option){
      if(fDebug){
        printf("    returning 1\n");
      }
      return 1;
    }
  }
  if(fDebug){
    printf("    returning 0\n");
  }
  return 0;
}


*/
