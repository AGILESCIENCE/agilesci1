////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_fitpsfarray
//       Release: BUILD 0.1 -  25/May/2010
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
//       25/May/2010
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
#include "pil.h"
#include "fitsio.h"
#include "CalibUtils.h"
#include "MathUtils.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TF1.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TFile.h"

using namespace std;

	Double_t fitking(Double_t *x, Double_t *params) {
//		TFormula f("kingfunction", "(1. - 1./[1]) * std::pow(1. +  ((x/[0])**2.)/(2.0*[1]), -[1])");
		Double_t f1 = (1. - 1./params[2]) * pow(1. +  (x[0]*x[0]/params[1]/params[1])/(2.0 * params[2]), -params[2]);
		Double_t f2 = (1. - 1./params[5]) * pow(1. +  (x[0]*x[0]/params[4]/params[4])/(2.0 * params[5]), -params[5]);
		return params[0] * f1 + params[3] * f2;
	}

void AG_fitpsfarray(char *  outfilename, char *  psdfilename) {
	
	TString basename(outfilename);

// read PSD matrix
	/* float ******psfgrid;
	float **psfrho;
	float **psfpsi;
	float **psftheta;
	float **psfphi;
	float **psfenergy;
	long naxes[5];
		psfgrid = new  float ***** [1];
		psfrho = new  float * [1];
		psfpsi = new  float * [1];		
		psftheta = new  float * [1];
		psfphi = new  float * [1];
		psfenergy = new  float * [1];
		alikeReadPsf(psdfilename, psfgrid, psfrho, psfpsi, psftheta, psfphi, psfenergy, naxes);
	 */
	PsfGrid psf(psdfilename);
	
	long naxes[5];
	naxes[0] = psf.Rhos().Size();
	naxes[1] = psf.Psis().Size();
	naxes[2] = psf.Energies().Size();
	naxes[3] = psf.Thetas().Size();
	naxes[4] = psf.Phis().Size();
	float * psfrho = new float[naxes[0]];
	
	float * psfene2 = new float[naxes[2]+1];
	for (int i=0;i<naxes[2];i++) psfene2[i]=psf.Energies()[i];
	psfene2[naxes[2]]=50000;
	
	float drho = psf.Rhos()[1]-psf.Rhos()[0];
	for (int rho=0; rho<naxes[0] ; rho++) psfrho[rho] = psf.Rhos()[rho]-0.5*drho ;
	TF1 * Tfitking = new TF1("Tfitking", fitking, psfrho[0], psfrho[naxes[0]-1]+drho, 6);

	ofstream ofout((basename+".txt").Data());

// prepare histograms
	TH2F norm[2];
	norm[0] = TH2F("Norm0","Norm 0;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	norm[1] = TH2F("Norm45","Norm 45;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	TH2F ang[2];
	ang[0] = TH2F("Ang0","Angular Scale 0;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	ang[1] = TH2F("Ang45","Angular Scale 45;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	TH2F gamma[2];
	gamma[0] = TH2F("Gamma0","Gamma 0;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	gamma[1] = TH2F("Gamma45","Gamma 45;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	TH2F chi2[2];
	chi2[0] = TH2F("Chi^2 0","Chi^2 0;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	chi2[1] = TH2F("Chi^2 45","Chi^2 45;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	TH2I ndf[2];
	ndf[0] = TH2I("NDF 0","NDF 0;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());
	ndf[1] = TH2I("NDF 45","NDF 45;E (MeV);#theta (deg)",naxes[2],psfene2,naxes[3]-1,psf.Thetas().Buffer());

	delete [] psfene2;
	gStyle->SetOptStat("");
	gStyle->SetOptFit(0000);
//	int badgammas = 0;

	for (int energy=0; energy<naxes[2]; energy++)
		for (int theta=0; theta<naxes[3]; theta++)
			for (int phi=0; phi<2; phi++) {
				TCanvas * c0 = new TCanvas("","");
				TString outname(basename);
				outname += "_";
				outname += int(psf.Energies()[energy]);
				outname += "_";
				outname += int(psf.Thetas()[theta]);
				outname += "_";
				outname += int(psf.Phis()[phi]);
				cout << outname << endl;

				TH1F psfhist("psf",(outname+";#theta (deg)").Data(),naxes[0]-1,psfrho);
				c0->SetBottomMargin(0.2);
				c0->SetFillColor(kWhite);
//				c0->SetLineColor(kBlack);
				psfhist.GetXaxis()->SetLabelSize(0.06);
				psfhist.GetYaxis()->SetLabelSize(0.06);
				psfhist.GetXaxis()->SetTitleSize(0.06);
				psfhist.GetYaxis()->SetTitleSize(0.06);
				psfhist.Set(naxes[0], psf.Values().Buffer()+psf.Values().AbsIndex(phi,theta,energy,0,0));
				float scale = 100.0/psfhist.GetBinContent(1);
				psfhist.Scale(scale);
				psfhist.GetXaxis()->SetRange(0,0);
				psfhist.Draw("HIST");
				Tfitking->SetParameter(0,100.0);
				if (energy >= 5) 
					Tfitking->SetParameter(1,0.1);
				else
					Tfitking->SetParameter(1,1.0);					
				Tfitking->SetParameter(2,1.5);
				Tfitking->SetParLimits(0,0.0,1000.0);
				Tfitking->SetParLimits(1,0.0,12.0);
				Tfitking->SetParLimits(2,0.0,12.0);
				Tfitking->FixParameter(3,0.0);
				Tfitking->FixParameter(4,0.1);
				Tfitking->FixParameter(5,0.1);
//				psfhist.Fit("Tfitking","WW","E0");
				Tfitking->ReleaseParameter(0);
				Tfitking->ReleaseParameter(1);
				Tfitking->ReleaseParameter(2);
//				psfhist.Fit("Tfitking","EMWW");
//				Tfitking->ReleaseParameter(3);
//				Tfitking->ReleaseParameter(4);
//				Tfitking->ReleaseParameter(5);
//				Tfitking->SetParameter(3,1.0);
//				psfhist.Fit("Tfitking","EMWW");

				ofout << psf.Energies()[energy] << " " << psf.Thetas()[theta] << " " << psf.Phis()[phi] << " " 
					<< Tfitking->GetParameter(0)/scale << " " << Tfitking->GetParameter(1) << " " << Tfitking->GetParameter(2) << " " <<  Tfitking->GetChisquare() << " " << endl;
					
				int maxnonzero = naxes[0];
				for ( ; psfhist.GetBinContent(maxnonzero) <= 5 ; maxnonzero--)
				psfhist.GetXaxis()->SetRange(1, maxnonzero);
				cout << maxnonzero << endl;
				outname += ".eps";
				c0->Print(outname.Data());
				delete c0;

				norm[phi].SetBinContent(energy+1,theta+1,Tfitking->GetParameter(0)/scale);
				norm[phi].SetBinError(energy+1,theta+1,Tfitking->GetParError(0)/scale);
				ang[phi].SetBinContent(energy+1,theta+1,Tfitking->GetParameter(1));
				ang[phi].SetBinError(energy+1,theta+1,Tfitking->GetParError(1));
				float p2=Tfitking->GetParameter(2);
//				if (p2 < 1e4) {
					gamma[phi].SetBinContent(energy+1,theta+1,p2);
					gamma[phi].SetBinError(energy+1,theta+1,Tfitking->GetParError(2));
//				} else {
//					gamma[phi].SetBinContent(energy+1,theta+1,1.0);
//					gamma[phi].SetBinError(energy+1,theta+1,Tfitking->GetParError(2));
//					badgammas++;
//				}
				chi2[phi].SetBinContent(energy+1,theta+1,Tfitking->GetChisquare());
				ndf[phi].SetBinContent(energy+1,theta+1,Tfitking->GetNDF());
	}
	ofout.close();
	delete [] psfrho;
	
	int tfields = 11;
	char ** ttype = new char*[tfields];
	char ** tform = new char*[tfields];
	char ** tunit = new char*[tfields];
	for (int i=0; i<tfields; i++) {
		ttype[i] = new char[FLEN_FILENAME];
		tform[i] = new char[FLEN_FILENAME];
		tunit[i] = new char[FLEN_FILENAME];
	}
	strcpy(ttype[0],"ENERGY");
	strcpy(tform[0],"1E");
	strcpy(tunit[0],"MeV");
	strcpy(ttype[1],"THETA");
	strcpy(tform[1],"1E");
	strcpy(tunit[1],"degrees");
	strcpy(ttype[2],"PHI");
	strcpy(tform[2],"1E");
	strcpy(tunit[2],"degrees");
	strcpy(ttype[3],"NORM");
	strcpy(tform[3],"1E");
	strcpy(tunit[3],"degrees^-1");
	strcpy(ttype[4],"NORM_ERR");
	strcpy(tform[4],"1E");
	strcpy(tunit[4],"degrees^-1");
	strcpy(ttype[5],"ANG");
	strcpy(tform[5],"1E");
	strcpy(tunit[5],"degrees");
	strcpy(ttype[6],"ANG_ERR");
	strcpy(tform[6],"1E");
	strcpy(tunit[6],"degrees");
	strcpy(ttype[7],"GAMMA");
	strcpy(tform[7],"1E");
	strcpy(tunit[7],"");
	strcpy(ttype[8],"GAMMA_ERR");
	strcpy(tform[8],"1E");
	strcpy(tunit[8],"");
	strcpy(ttype[9],"CHI2");
	strcpy(tform[9],"1E");
	strcpy(tunit[9],"");
	strcpy(ttype[10],"NDF");
	strcpy(tform[10],"1J");
	strcpy(tunit[10],"");
	
	fitsfile * fp;
	int status=0;
	fits_create_file(&fp, (basename+"_table.fits.gz").Data(), &status);
	if (status==0) fits_create_tbl(fp, BINARY_TBL, 0, tfields, ttype, tform, tunit, "PSDPARAMS", &status);
	int rownum = 0;
	if (status==0) for (int energy=0; energy<naxes[2]; energy++)
		for (int theta=0; theta<naxes[3]; theta++)
			for (int phi=0; phi<2; phi++) {
				float value = psf.Energies()[energy];
				if (status==0) fits_write_col(fp, TFLOAT, 1, ++rownum, 1, 1, &value, &status);
				value = psf.Thetas()[theta];
				if (status==0) fits_write_col(fp, TFLOAT, 2, rownum, 1, 1, &value, &status);
				value = psf.Phis()[phi];
				if (status==0) fits_write_col(fp, TFLOAT, 3, rownum, 1, 1, &value, &status);
				value = norm[phi].GetBinContent(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 4, rownum, 1, 1, &value, &status);
				value = norm[phi].GetBinError(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 5, rownum, 1, 1, &value, &status);
				value = ang[phi].GetBinContent(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 6, rownum, 1, 1, &value, &status);
				value = ang[phi].GetBinError(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 7, rownum, 1, 1, &value, &status);
				value = gamma[phi].GetBinContent(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 8, rownum, 1, 1, &value, &status);
				value = gamma[phi].GetBinError(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 9, rownum, 1, 1, &value, &status);
				value = chi2[phi].GetBinContent(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TFLOAT, 10, rownum, 1, 1, &value, &status);
				int value2 = ndf[phi].GetBinContent(energy+1,theta+1);
				if (status==0) fits_write_col(fp, TLONG, 11, rownum, 1, 1, &value2, &status);
	}
	if (status==0) fits_close_file(fp, &status);
	
	TCanvas * c0 = new TCanvas("","");
	gPad->SetLogx(1);
	TFile f(basename+"_histos.root","recreate");
	for (Long_t phi=0;phi<2;phi++){
		gPad->SetLogz(1);
		norm[phi].Draw("SURF1");
		c0->Print((basename+phi+"_norm.eps").Data());
		norm[phi].Write();
		ang[phi].Draw("SURF1");
		c0->Print((basename+phi+"_ang.eps").Data());
		ang[phi].Write();
		gPad->SetLogz(0);
		gamma[phi].SetMaximum(4.5);
		gamma[phi].Draw("SURF1");
		c0->Print((basename+phi+"_gamma.eps").Data());
		gamma[phi].Write();
	}
	delete c0;
//	cerr << badgammas << " bad gammas" << endl;
}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char * outfilename = new char[FLEN_FILENAME];
	char * psdfilename  = new char[FLEN_FILENAME];
	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("outfile", outfilename);
	status = PILGetString("psdfile", psdfilename);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_fitpsfarray.cpp v.0.1 -25/05/10 - A.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Output file name : " <<  outfilename << endl;
	cout << "PSD file name : "<< psdfilename << endl;
	

	cout << "AG_fitpsfarray...............................starting"<< endl;		
	AG_fitpsfarray(outfilename, psdfilename);
	cout << "AG_fitpsfarray............................... exiting"<< endl;		
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_writepsf........... exiting #################\n");
		printf("#################################################################\n\n\n");					
	
	delete[] outfilename;
	delete[] psdfilename;

	return status;
}

