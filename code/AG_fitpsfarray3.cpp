////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_fitpsfarray3
//       Release: BUILD 0.1 -  10/Jun/2010
//       Contributors: 
//       Author: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       10/Jun/2010
//                      First release: V0.1
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "pil.h"
#include "fitsio.h"
#include "MathUtils.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TFile.h"

using namespace std;

const Double_t drho = 0.1;
const int nrhos = 200;
const Double_t maxrho = drho * nrhos;
const Double_t normscale = 1e5;
const Double_t psfcoeff = M_PI * M_PI / 90.0 * drho * normscale;

Double_t fitking(Double_t *x, Double_t *params) {
//		TFormula f("kingfunction", "(1. - 1./[1]) * std::pow(1. +  ((x/[0])**2.)/(2.0*[1]), -[1])");
	Double_t f1 = (1. - 1./params[2]) * pow(1. +  (x[0]*x[0]/params[1]/params[1])/(2.0 * params[2]), -params[2]);
	Double_t f2 = (1. - 1./params[5]) * pow(1. +  (x[0]*x[0]/params[4]/params[4])/(2.0 * params[5]), -params[5]);
	return psfcoeff * sin(x[0] * DEG2RAD) * (params[0] * f1 + params[3] * f2);
}

class FitResults{
public:
	Double_t norm, norm_err, gamma, gamma_err, ang, ang_err, chi2, ndf;

	FitResults():norm(0),norm_err(0), gamma(0), gamma_err(0), ang(0), ang_err(0), chi2(0), ndf(0){}

	void Set(Double_t _norm, Double_t _norm_err, Double_t _gamma, Double_t _gamma_err, Double_t _ang, Double_t _ang_err, Double_t _chi2, Double_t _ndf){
		norm = _norm;
		norm_err = _norm_err;
		gamma = _gamma;
		gamma_err = _gamma_err;
		ang = _ang;
		ang_err = _ang_err;
		chi2 = _chi2;
		ndf = _ndf;
	}

	void FillFromFunction(TF1& infunc){
		norm = normscale * infunc.GetParameter(0);
		norm_err = normscale * infunc.GetParError(0);
		gamma = infunc.GetParameter(2);
		gamma_err = infunc.GetParError(2);
		ang = infunc.GetParameter(1);
		ang_err = infunc.GetParError(1);
		chi2 = infunc.GetChisquare();
		ndf = infunc.GetNDF();
	}

};

ostream& operator<<(ostream& _stream, FitResults& f){
	_stream << " " << f.norm << " " << f.norm_err << " " << f.gamma << " " << f.gamma_err << " " << f.ang << " " << f.ang_err << " " << f.chi2 << " " << f.ndf;
	return _stream;
}
	

void AG_fitpsfarray(char *  dataprefix, char * theta, char * phi) {
	
	TCanvas * c0 = new TCanvas("","");
	c0->SetFillColor(kWhite);
	c0->SetBottomMargin(0.2);
	int itheta = atoi(theta);
	int iphi = atoi(phi);
	TString dpre(TString(dataprefix)+"_" + theta + "_" + phi);
	TString labelpre(TString("SIM")+"_" + theta + "_" + phi);
	
	gStyle->SetOptStat("ne");
	gStyle->SetOptFit(111);

	const int nfilters = 3;
	const TString filter[] = {"_FM","_F4","_FT3ab"};
	const int nenergies = 16;
	const float psfenergy[] = {10.,35,50,71,100,141,200,283,400,632,1000,1732,3000,5477,10000,20000,100000};
	const int nevtypes = 3;
	const char evtype[] = "GLS";

	TF1 * Tfitking = new TF1("Tfitking", fitking, 0.0, maxrho, 6);
	Tfitking->SetParName(0,"B");
	Tfitking->SetParName(1,"#delta");
	Tfitking->SetParName(2,"#gamma");
	
//read event energies
	vector<int> evcanals;
	TString datfilename(dpre+".dat");
	cout << "Reading " << datfilename << endl;
	ifstream datfile(datfilename);
	long evnum = 0;
	long ensum[nenergies];
	cout << "Minimum energy = " << psfenergy[0] << endl;
	cout << "Maximum energy = " << psfenergy[nenergies] << endl;
	for (int i=0; i<nenergies ; i++) ensum[i]=0;
	for ( ; !datfile.eof() ; ) {
		long dummy;
		float evenergy;
		string dumstring;
		datfile >> dummy >> evenergy;
		if (!datfile.eof()) {
			getline(datfile, dumstring);
			evenergy *= 1000;
			if (evenergy < psfenergy[0])
				evcanals.push_back(-1);
			else if (evenergy < psfenergy[nenergies]) {
				int canal=0;
				for ( ; evenergy >= psfenergy[canal] ; canal++) {}
				canal--;
				evcanals.push_back(canal);
				ensum[canal]++;
			} else
				evcanals.push_back(nenergies);
			evnum++;
		}
	}
	datfile.close();
	cout << evnum << " events" << endl;
	for (int i=0; i<nenergies; i++) cout << ensum[i] << " events in channel " << i << endl;
	
	for (int f=0;f<nfilters;f++){
		FitResults psfresults[nenergies*nevtypes];
// open flg file
		int status=0;
		fitsfile * flg;
		TString flgfilename(dpre + filter[f] + ".flg");
		cout << "FLG file name = " << flgfilename << endl;
		if ( fits_open_table(&flg, flgfilename.Data(), READONLY, &status) !=0 ) {
			cerr << "Could not open " << flgfilename << endl;
		} else {
// read flg file
			int thetacol, phicol, evcol;
			fits_get_colnum(flg, CASEINSEN, "THETA", &thetacol, &status);
			fits_get_colnum(flg, CASEINSEN, "PHI", &phicol, &status);
			fits_get_colnum(flg, CASEINSEN, "EVSTATUS", &evcol, &status);
			long testnrows;
			fits_get_num_rows(flg, &testnrows, &status);
			if (evnum != testnrows) cerr << "Warning: " << evnum << " events in " << datfilename << " do not match " << testnrows << " events in " << flgfilename << endl;
			float * evtheta = new float[evnum];
			float * evphi = new float[evnum];
			char * evevtype = new char[evnum];
			fits_read_col(flg, TFLOAT, thetacol, 1, 1, testnrows, NULL, evtheta, NULL, &status);
			fits_read_col(flg, TFLOAT, phicol, 1, 1, testnrows, NULL, evphi, NULL, &status);
			fits_read_col(flg, TBYTE, evcol, 1, 1, testnrows, NULL, evevtype, NULL, &status);
			fits_close_file(flg, &status);
			
// Initialize histograms
			TH1I psfhist[nenergies*nevtypes];
			for (int en=0; en<nenergies ; en++) for (int ev=0; ev<nevtypes ; ev++) {
				long index = ev*nenergies+en;
//				psfhist[index] =  TH1I(TString(evtype[ev])+Long_t(psfenergy[en]),(labelpre+filter[f]+evtype[ev]+Long_t(psfenergy[en])+";#theta (deg)").Data(),nrhos,0.0,maxrho); 
				psfhist[index] =  TH1I(TString("King Function"),(labelpre+filter[f]+evtype[ev]+Long_t(psfenergy[en])+";#theta (deg)").Data(),nrhos,0.0,maxrho); 
				psfhist[index].GetXaxis()->SetLabelSize(0.06);
				psfhist[index].GetYaxis()->SetLabelSize(0.06);
				psfhist[index].GetXaxis()->SetTitleSize(0.06);
				psfhist[index].GetYaxis()->SetTitleSize(0.06);
				psfhist[index].Sumw2();
			}
// Fill histograms
			for (long row = 0 ; row < evnum ; row++ ) {
				int ev = 0;
				for ( ;  ev < 3 && evtype[ev] != evevtype[row]; ev++) {}
				if (ev < 3 && evcanals[row] >=0 && evcanals[row]<nenergies)
					if (f != 2 || evtheta[row] != (-1.0) || evphi[row] != (-1.0) || ev != 2)
						psfhist[ev*nenergies+evcanals[row]].Fill(SphDistDeg(evphi[row], 90.0-evtheta[row], iphi, 90.0-itheta));
			}
// Fit histograms
			 for (int ev=0; ev<nevtypes ; ev++) {
				TString psffilename(dpre+filter[f]+evtype[ev]+".psf3");
				ofstream ofout(psffilename.Data());
				cout << "Writing " << psffilename << endl;
				for (int en=0; en<nenergies ; en++){
					long index = ev*nenergies+en;
					cout << psfenergy[en] << " MeV = channel " << en << endl;
					Double_t psfsum = psfhist[index].GetSum();
					cout << psfsum << " events" << endl;
					if (psfsum >= 10) {
						psfhist[index].GetXaxis()->SetRange(0,0);
						Tfitking->SetParameter(0, 1e-3 * psfsum);
						if (en >= 5) 
							Tfitking->SetParameter(1,0.1);
						else
							Tfitking->SetParameter(1,1.0);					
						Tfitking->SetParameter(2,1.5);
						Tfitking->SetParLimits(0,0.0, 0.1 * psfsum);
						Tfitking->SetParLimits(1,0.0,12.0);
						Tfitking->SetParLimits(2,0.0,12.0);
						Tfitking->FixParameter(3,0.0);
						Tfitking->FixParameter(4,0.1);
						Tfitking->FixParameter(5,0.1);
						psfhist[index].Fit("Tfitking","LL","E0");
						Tfitking->ReleaseParameter(0);
						Tfitking->ReleaseParameter(1);
						Tfitking->ReleaseParameter(2);
						psfhist[index].Fit("Tfitking","LLEM");
						psfresults[index].FillFromFunction(*Tfitking);
						ofout << psfenergy[en] << " " << psfresults[index] << endl;
							
						int xlength = nrhos;
						Int_t histmin = 0.05 * psfhist[index].GetMaximum();
						for ( ; psfhist[index].GetBinContent(xlength) <= histmin ; xlength--)
						psfhist[index].GetXaxis()->SetRange(1, xlength);
						c0->Print((dpre+filter[f]+evtype[ev]+Long_t(psfenergy[en])+".eps").Data());
					}
				}
				ofout.close();
			}
		}
	}
		
	delete c0;
	delete Tfitking;
//	cerr << badgammas << " bad gammas" << endl;
}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char * dataprefix = new char[FLEN_FILENAME];
	char * theta = new char[2];
	char * phi = new char[2];
	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("dataprefix", dataprefix);
	status = PILGetString("theta", theta);
	status = PILGetString("phi", phi);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_fitpsfarray3.cpp v.0.1 -25/05/10 - A.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "data file prefix : "<< dataprefix << endl;
	cout << "theta : "<< theta << endl;
	cout << "phi : "<< phi << endl;
	

	cout << "AG_fitpsfarray...............................starting"<< endl;		
	AG_fitpsfarray(dataprefix, theta, phi);
	cout << "AG_fitpsfarray............................... exiting"<< endl;		
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_fitpsfarray........... exiting #################\n");
		printf("#################################################################\n\n\n");					
	
	delete[] theta;
	delete[] phi;
	delete[] dataprefix;

	return status;
}

