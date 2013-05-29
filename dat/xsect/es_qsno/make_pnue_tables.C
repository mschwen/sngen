//
// make_pnue_tables.C
//
// make   qphysics/parameters/pnue_tot.dat
//        qphysics/parameters/pnue_dif.dat
//
//
{
    gROOT->Reset();

    //-------------------------------------------------------------------
    // Set table parameters here 
    //-------------------------------------------------------------------
    Int_t    MakeTotal  = 0;       //  =1 make,  =0 not make
    Double_t EnuStepTot = 0.1;
    Int_t    NDataTot   = 1700;    // EnuMax = EnuStep * NData

    Int_t    MakeDiff   = 1;       //  =1 make,  =0 not make
    Double_t EnuStepDif = 0.1;
    Int_t    NDataDif   = 1700;    // EnuMax = EnuStep * NData

    Double_t me      = 0.510998918; // MeV; updated to PDG2004 value by JH

    //-------------------------------------------------------------------
    // open output files
    //-------------------------------------------------------------------
    if (MakeTotal == 1) {
	FILE *fp_tot = fopen("pnue_tot.dat", "w"); 
	fprintf(fp_tot, "EnuStepTot = %f\n", EnuStepTot);
	fprintf(fp_tot, "NDataTot   = %d\n", NDataTot);
    }
    if (MakeDiff == 1) {
	FILE *fp_dif = fopen("pnue_dif.dat", "w"); 
	fprintf(fp_dif, "EnuStepDif = %f\n", EnuStepDif);
	fprintf(fp_dif, "NDataDif   = %d\n", NDataDif);
    }

    //-------------------------------------------------------------------
    // prepare PnuE function
    //-------------------------------------------------------------------
    PNuE *pnue;
    pnue = new PNuE(0);               // don't initialize
    pnue->SetRadiativeCorrection(3);  // radiative correction w/o table
    pnue->SetVerbosity(0);            // minimum messages to save time

    //-------------------------------------------------------------------
    // initialize clock to monitor running time (optional)
    //-------------------------------------------------------------------
#include <time.h>
    clock_t start = clock();

    //-------------------------------------------------------------------
    // main loop (total cross section)
    //-------------------------------------------------------------------

    Double_t Enu, Te, Temax, sigma_nue, sigma_numu;

    if (MakeTotal == 1) {
	for (Int_t i = 1; i <= NDataTot; i++) {
	    Enu =  EnuStepTot * (double) i;
	    pnue->SetReaction("nue");
	    sigma_nue = pnue->Sigma(Enu);
	    pnue->SetReaction("numu");
	    sigma_numu = pnue->Sigma(Enu);
	    // output
	    fprintf(fp_tot, "  %6.3lf  %20.12le  %20.12le\n", 
		    Enu, sigma_nue, sigma_numu);
	    // monitor
	    if (i%10 == 0) {
		printf("  %6.3lf  %20.12le  %20.12le\n", 
		       Enu, sigma_nue, sigma_numu);
	    }
	}
	fclose(fp_tot);
    }

    //-------------------------------------------------------------------
    // check & reset clock (optional)
    //-------------------------------------------------------------------
    double tsec = (double)(clock()-start) / 1000000.;
    printf("pnue_tot.dat: rough time = %f sec  loop=%g  ave=%g sec\n", 
	   tsec, NDataTot,  
	   tsec / (double)NDataTot); 
    start = clock();

    //-------------------------------------------------------------------
    // main loop (differential cross section)
    //-------------------------------------------------------------------
    if (MakeDiff == 1) {
	for (Int_t i = 1; i <= NDataDif; i++) {
	    Enu =  EnuStepDif * (double) i;
	    Temax = Enu - 1.0/(2.0/me + 1.0/Enu);
	    for (Int_t j = 1; j < NDataDif; j++) {
		Te =  EnuStepDif * (double) j;
		if (Te <= Temax) {
		    pnue->SetReaction("nue");
		    sigma_nue = pnue->dSigmadT(Enu, Te);
		    pnue->SetReaction("numu");
		    sigma_numu = pnue->dSigmadT(Enu, Te);
		    // output
		    fprintf(fp_dif, "  %6.3lf  %6.3lf  %20.12le  %20.12le\n", 
			    Enu, Te, sigma_nue, sigma_numu);
		}
	    }
	    // monitor
	    if (Te <= Temax && i%10 == 0) 
		printf("  %6.3lf  %6.3lf  %20.12le  %20.12le\n", 
		       Enu, Te, sigma_nue, sigma_numu);
	}
	fclose(fp_dif);
    }

    //-------------------------------------------------------------------
    // check clock (optional)
    //-------------------------------------------------------------------
    tsec = (double)(clock()-start) / 1000000.;
    printf("pnue_dif.dat: rough time = %f sec  loop=%g  ave=%g sec\n", 
	   tsec, NDataDif,  
	   tsec / (double)NDataDif); 

}
