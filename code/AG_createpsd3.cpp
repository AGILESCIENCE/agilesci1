////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_createpsd3
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
//       20/Jun/2010
//                      First release: V0.1
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include "pil.h"
#include "fitsio.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

const double drho = 0.1;
const int numrhos = 200;
const double maxrho = drho * numrhos;
const double normscale = 1e5;
const double dtor = M_PI / 180.0;
const double psfcoeff = M_PI * M_PI / 90.0 * drho;
const char * thetas[] = {"01","30","35","45","50","60"};
const double fthetas[] = {1.0, 30.0, 35.0, 45.0, 50.0, 60.0};
const char * phis[] = {"00","45"};
const int numthetas = 6;
const int numphis = 2;
const float dtheta = 5.0;
const float dphi = 45.0;
const float numdthetas = 19;
const float numdphis = 8;

double fitking(double *x, double *params) {
//		TFormula f("kingfunction", "(1. - 1./[1]) * std::pow(1. +  ((x/[0])**2.)/(2.0*[1]), -[1])");
	double f1 = (1. - 1./params[2]) * pow(1. +  (x[0]*x[0]/params[1]/params[1])/(2.0 * params[2]), -params[2]);
	return params[0] * f1;
}

class FitResults{
public:
	double norm, norm_err, gamma, gamma_err, ang, ang_err, chi2, ndf;

	FitResults():norm(0),norm_err(0), gamma(0), gamma_err(0), ang(0), ang_err(0), chi2(0), ndf(0){}
	~FitResults(){}
	
	void Set(double _norm, double _norm_err, double _gamma, double _gamma_err, double _ang, double _ang_err, double _chi2, double _ndf){
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

class FitResultList{
public:
	vector<float> energies;
	vector<FitResults> fitresults;

	FitResultList(){}
	~FitResultList(){}
	
	bool Read(const char * filename){
		ifstream infile(filename);
		cout << filename << endl;
		bool fail = infile.fail();
		for ( ; !infile.eof() ; ) {
			float energy;
			FitResults fitresult;
			infile >> energy;
			if (!infile.eof()) {
				infile >> fitresult.norm >> fitresult.norm_err 
				>> fitresult.gamma >> fitresult.gamma_err 
				>> fitresult.ang >> fitresult.ang_err
				>> fitresult.chi2 >> fitresult.ndf ;
				fail = fail || infile.fail();
				energies.push_back(energy);
				fitresults.push_back(fitresult);
//				cout << energies.back() << " " << fitresults.back().gamma << " " << fitresults.back().ang << endl;
			}
		}
		infile.close();
		if (fail) cerr << "fail" << endl;
		return fail;
	}
	
	FitResultList(const char * filename){ Read(filename); }
};

class ParamGrid{
public:

	FitResultList fitlists[numthetas * numphis];
	TString _filterev;
	double rhos[numrhos];
	
	ParamGrid():_filterev("_FMG"){}
	~ParamGrid(){}
	
	inline long ind(int t, int p){return t*numphis+p;}

	bool Read(const char * prefix, const char * filterev = "_FMG") {
		_filterev = filterev;
		bool fail = false;
		for (int i=0; i<numrhos; i++) rhos[i] = drho * (0.5 + i);
		for (int t=0; t<numthetas ; t++) for (int p=0; p<numphis; p++) {
			TString filename = ((TString(prefix) + "_") + thetas[t]) + "_" ;
			if (t==0) filename += phis[0]; else filename += phis[p];
			filename += _filterev + ".psf3";
			fail = fail || fitlists[ind(t,p)].Read(filename.Data());
		}
		return fail;
	}
	
	ParamGrid(const char * prefix, const char * filterev = "_FMG"){Read(prefix, filterev);}

	
//	inline FitResults& operator()(int e, int t, int p){return fitlists[ind(t,p)].fitresults[e];}

	void findep(float energy, float phi, int& e, int&p){
		e = 0;
		int maxe = fitlists[0].energies.size() - 1;
		for ( ; fitlists[0].energies[e] < energy ; e++){};
		if (e > 0) e--;
		if (e > maxe) e = maxe;
		p = 0;
		phi = fmod(phi + 720.0, 360.0);
		for ( ; p*dphi < phi ; p++){}
		p %= 2;
	}

	int findt(float theta){
		int t = 0;
		for ( ; t < numthetas && fthetas[t] <= theta ; t++){}
		if (t > 0) t--;
		if (t > numthetas-1) t = numthetas-1;
		return t;
	}
	
	double angle(int e, int p, float theta){
		double angles[numthetas];
		for (int t=0; t<numthetas; t++)
			angles[t] = fitlists[ind(t,p)].fitresults[e].ang;
		if (theta > fthetas[numthetas-1])
			return angles[numthetas-1];
		else {
			TGraph anggraph(numthetas, fthetas, angles);
			return anggraph.Eval(theta);
		}
	}

	double gamma(int e, int p, float theta){
		double gammas[numthetas];
		for (int t=0; t<numthetas; t++)
			gammas[t] = fitlists[ind(t,p)].fitresults[e].gamma;
		if (theta > fthetas[numthetas-1])
			return gammas[numthetas-1];
		else {
			TGraph gammagraph(numthetas, fthetas, gammas);
			return gammagraph.Eval(theta);
		}
	}

	double norm(double_t gamma, Double_t ang){
		double params[3];
		params[0] = 1.0;
		params[1] = ang;
		params[2] = gamma;
		double psf[numrhos];
		double psfsum=0.0;
		for (int i=0; i<numrhos; i++){
			psf[i] = fitking(&rhos[i], params);
			psfsum += psf[i] * sin(rhos[i] * dtor);
		}
		return 1.0 / (psfcoeff * psfsum);
	}
	
	int Write(const char * outfilename){
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
		fits_create_file(&fp, outfilename, &status);
		if (status==0) 
			fits_create_tbl(fp, BINARY_TBL, 0, tfields, ttype, tform, tunit, "PSDPARAMS", &status);
		else
			cerr << "File not created : status = " << status << endl;
		long rownum = 0;
		cout << fitlists[0].energies.size() << " " ;
		if (status==0) for (unsigned int e=0; e < fitlists[0].energies.size(); e++) {
			float energy = fitlists[0].energies[e];
			for (float tt=0; tt<numdthetas; tt++) {
				float theta = 0.0 + tt * dtheta;
				int t = findt(theta);
				for (int pp=0; pp<numphis; pp++) {
					float phi = 0.0 + pp * dphi;
					int p = pp % 2;
					cout << energy << " " << theta << " " << phi << " ";
					cout << e << " " << t << " " << p << " ";
					if (status==0) fits_write_col(fp, TFLOAT, 1, ++rownum, 1, 1, &energy, &status);
					if (status==0) fits_write_col(fp, TFLOAT, 2, rownum, 1, 1, &theta, &status);
					if (status==0) fits_write_col(fp, TFLOAT, 3, rownum, 1, 1, &phi, &status);
					float angval = angle(e, p, theta);
					if (status==0) fits_write_col(fp, TFLOAT, 6, rownum, 1, 1, &angval, &status);
					cout << angval << endl;
					float value = angval * fitlists[ind(t,p)].fitresults[e].ang_err / fitlists[ind(t,p)].fitresults[e].ang;
					if (status==0) fits_write_col(fp, TFLOAT, 7, rownum, 1, 1, &value, &status);
					float gammaval = gamma(e, p, theta);
					if (status==0) fits_write_col(fp, TFLOAT, 8, rownum, 1, 1, &gammaval, &status);
					value = gammaval * fitlists[ind(t,p)].fitresults[e].gamma_err / fitlists[ind(t,p)].fitresults[e].gamma;
					if (status==0) fits_write_col(fp, TFLOAT, 9, rownum, 1, 1, &value, &status);
					float normval = norm(gammaval, angval);
					if (status==0) fits_write_col(fp, TFLOAT, 4, rownum, 1, 1, &normval, &status);
					value = normval * fitlists[ind(t,p)].fitresults[e].norm_err / fitlists[ind(t,p)].fitresults[e].norm;
					if (status==0) fits_write_col(fp, TFLOAT, 5, rownum, 1, 1, &value, &status);
					value = fitlists[ind(t,p)].fitresults[e].chi2;
					if (status==0) fits_write_col(fp, TFLOAT, 10, rownum, 1, 1, &value, &status);
					value = fitlists[ind(t,p)].fitresults[e].ndf;
					if (status==0) fits_write_col(fp, TLONG, 11, rownum, 1, 1, &value, &status);
				}
			}
		} else
			cerr << "Table not created: status = " << status << endl;
		if (status==0) fits_close_file(fp, &status);
		for (int i=0; i<tfields; i++) {
			delete [] ttype[i];
			delete [] tform[i];
			delete [] tunit[i];
		}
		delete [] ttype;
		delete [] tform;
		delete [] tunit;
		return status;
	}
};


int AG_createpsd3(const char *  dataprefix, const char * outprefix) {
	
	TString opre(outprefix);
	int status=0;
	const int nfilters = 3;
	const TString filter[] = {"_FM","_F4","_FT3ab"};
	const int nevtypes = 3;
	const char evtype[] = "GLS";
	
	for (int ev=0; ev<nevtypes; ev++) for (int f=0;f<nfilters;f++) {
		TString filterev = filter[f] + evtype[ev];
		TString outfilename = opre + filterev + "_table.fits.gz";
		ParamGrid paramgrid(dataprefix, filterev.Data());
		cout << "Writing to " << outfilename << endl;
		status = paramgrid.Write(outfilename.Data());
	}
	return status;
}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char * dataprefix  = new char[FLEN_FILENAME];
	char * outprefix = new char[FLEN_FILENAME];
	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("dataprefix", dataprefix);
	status = PILGetString("outprefix", outprefix);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_createpsd3.cpp v.0.1 -20/06/10 - A.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Data filename prefix : "<< dataprefix << endl;
	cout << "Output filename prefix : " <<  outprefix << endl;
	

	cout << "AG_createpsd3...............................starting"<< endl;		
	status = AG_createpsd3(dataprefix, outprefix);
	if (status == 0)
		cout << "AG_createpsd3............................... exiting"<< endl;		
	else {
		cerr << "AG_createdpsd3: error = " << status << endl;
	}
	delete[] dataprefix;
	delete[] outprefix;

	printf("\n\n\n###################################################################\n");
	printf("#########  Task AG_createpsd3........... exiting #################\n");
	printf("#################################################################\n\n\n");					
		
	return status;
}

