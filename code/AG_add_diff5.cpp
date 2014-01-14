////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_add_diff
//       Release: V2.1 -  09/Oct/2010
//       Contributors: 
//       Author: Andrew Chen, Tomaso Contessi (IASF-Milano)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       27/Feb/2009
//                      First release: V1.0
//       		Author: Andrew Chen (IASF-Milano)
//       01/Oct/2010
//                      V2.0
//       		Author: Andrew Chen (IASF-Milano)
//       09/Oct/2010
//                      V2.1
//       		Author: Andrew Chen (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include "fitsio.h"
#include "pil.h"
#include <CalibUtils.h>

using namespace std;

/**
int AG_findindex(double emin, double * raeffenergy, int numaeffenergies) {
	int lowetrue = numaeffenergies - 1;
	if (raeffenergy[lowetrue] >= emin) {
	  for (lowetrue=0; raeffenergy[lowetrue] < emin; lowetrue++);
	  if (lowetrue > 0 && emin / raeffenergy[lowetrue-1] < raeffenergy[lowetrue] / emin)
		lowetrue--;
	}
	return lowetrue;
}
*/

int AG_add_diff(char * diffusefilelist, char * sarfile, char * edpfile, char *  outfile, double emin, double emax){
	
	EdpGrid edpgrid(edpfile);
	int status = 0;
	long pixel[2] = { 1, 1 };

	ifstream infile(diffusefilelist);
	double bigindex;
	int numdiffs;

	infile >> bigindex >> numdiffs;
	cout << "Index = " << bigindex << endl;

	/// AlikeAeffGridClass3 raeff(sarfile, edpfile, emin, emax, bigindex);
	AeffGridAverage raeff(sarfile, emin, emax, bigindex);
	raeff.LoadEdp(edpfile);

	/// double bigarea = raeff.Valavg(30,0);
	double bigarea = raeff.AvgVal(30, 0);

	cout << bigarea << endl;


	/**
	int numaeffenergies = raeff.energies().GetNoElements();
	double * raeffenergies = new double[numaeffenergies];
	for (int eind = 0 ; eind < numaeffenergies ; eind++)
	    raeffenergies[eind] = raeff.energies()[eind];
	int eminind = AG_findindex(emin, raeffenergies, numaeffenergies);
	int emaxind = AG_findindex(emax, raeffenergies, numaeffenergies);
	*/
	const VecF& raeffenergies = raeff.Energies();
	int numaeffenergies = raeffenergies.Size();
	int eminind = raeffenergies.GeomIndex(emin);
	int emaxind = raeffenergies.GeomIndex(emax);


	string diffusefilename;
	double * diffuseout;

	int bitpix   =  DOUBLE_IMG;
	int naxis;   
	long naxes[3] ;

	fitsfile * outFits;
	if ( fits_create_file(&outFits, outfile, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", outfile);
		return status;
	}

	for (int diffi = 0; diffi < numdiffs ; ++diffi){
	  infile >> diffusefilename;
	  fitsfile * diffuseFits;
	  if ( fits_open_file(&diffuseFits, diffusefilename.c_str(), READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", diffusefilename.c_str());
		return status;
		}
	  else {
	    fits_movabs_hdu(diffuseFits, 2, NULL, &status);	
	    fits_get_img_param(diffuseFits, 3, &bitpix, &naxis, naxes, &status);
	    cout << diffusefilename << ": " << naxis << " axes: (" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << ")" << endl;
	    if (diffi == 0) {
		fits_copy_file(diffuseFits+diffi, outFits, 1, 1, 1, &status);
		diffuseout = new double[naxes[0] * naxes[1]];
		for (long i=0 ; i < naxes[0] * naxes[1] ; i++)
		    diffuseout[i] = 0;
	    }
	
	    double * diffuse = new double[naxes[0] * naxes[1]];
	    if ( fits_read_pix(diffuseFits, TDOUBLE, pixel, naxes[0]*naxes[1], NULL, diffuse, NULL, &status) != 0) {
		printf("Error reading array from '%s'\n", diffusefilename.c_str());
		return status;
	    }
	    else {
		cout << "Diffuse array " << diffi << " for " << diffusefilename << " successfully created" << endl;
	
		double elow;
		double ehigh;
		double index;
		fits_read_key(diffuseFits,TDOUBLE,"E_MIN",&elow,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"E_MAX",&ehigh,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"INDEX",&index,NULL,&status);
		/**
		int elowind = AG_findindex(elow, raeffenergies, numaeffenergies);
		int ehighind = AG_findindex(ehigh, raeffenergies, numaeffenergies);
		*/
		int elowind = raeffenergies.GeomIndex(elow);
		int ehighind = raeffenergies.GeomIndex(ehigh);

		cout << "Elow = " << elow << ", Ehigh = " << ehigh << ", index = " << index << endl;

		/// AlikeAeffGridClass2 raeff2(sarfile, elow, ehigh, index);
		AeffGridAverage raeff2(sarfile, elow, ehigh, index);



		double specwttotal = 0;
		double edparr = 0;
		for (int etrue = elowind; etrue < ehighind; etrue++) {
		    cout << etrue << " " << raeffenergies[etrue] << endl;
		    double specwt;
		    if (etrue >= numaeffenergies)
			/// specwt = pow((double) raeffenergies[etrue] ,(double) (1.0 - index));
			specwt = pow(double(raeffenergies[etrue]), 1.0-index);
		    else
			/// specwt = pow((double) raeffenergies[etrue] ,(double) (1.0 - index)) - pow((double)raeffenergies[etrue+1], (double) (1.0 - index)) ;
			specwt = pow(double(raeffenergies[etrue]),1.0-index) - pow(double(raeffenergies[etrue+1]), 1.0 - index);


		    specwttotal += specwt;
		    for (int eobs = eminind;  eobs <= emaxind; eobs++) {
			edparr += specwt * edpgrid.Val(raeffenergies[etrue], raeffenergies[eobs], 30, 0);
		    cout << eobs << " " << raeffenergies[eobs] << " " << specwt << " " << edpgrid.Val(raeffenergies[etrue], raeffenergies[eobs], 30, 0) << endl;
		    }
		}
		cout << edparr / specwttotal << " " << raeff2.AvgVal(30,0) << endl;
		edparr *= raeff2.AvgVal(30,0) / specwttotal / bigarea;
		cout << edparr << endl;
		for (long i=0 ; i < naxes[0] * naxes[1] ; i++) {
		    diffuseout[i] += edparr * diffuse[i];
		}
		delete [] diffuse;
	    }
	  }
	  fits_close_file(diffuseFits, &status);
	}
	
	fits_movabs_hdu(outFits, 2, NULL, &status);
	fits_write_pix(outFits, TDOUBLE, pixel, naxes[0]*naxes[1], diffuseout, &status);	
	fits_update_key(outFits, TDOUBLE, "E_MIN", &emin, NULL, &status);
	fits_update_key(outFits, TDOUBLE, "E_MAX", &emax, NULL, &status);
	fits_update_key(outFits, TDOUBLE, "INDEX", &bigindex, NULL, &status);
	fits_close_file(outFits, &status);
	delete [] diffuseout;
	return status;
}

int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char diffusefilelist[FLEN_FILENAME];
	char sarfile[FLEN_FILENAME];
	char edpfile[FLEN_FILENAME];
	char outfile[FLEN_FILENAME];
	double emin = 100.0;
	double emax = 50000.0;
	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("diffusefilelist", diffusefilelist);	
	status = PILGetString("sarfile", sarfile);
	status = PILGetString("edpfile", edpfile);
	status = PILGetString("outfile", outfile);
	status = PILGetReal("emin", &emin);
	status = PILGetReal("emax", &emax);
	
	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_add_diff.cpp v.1 - 28/02/09 - A.C. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Enter diffuse file list = "<< diffusefilelist << endl;
	cout << "Enter sensitive area file name = "<< sarfile << endl;
	cout << "Enter energy dispersion file name = "<< edpfile << endl;
	cout << "Enter output file name = " <<outfile<< endl;
	cout << "Enter minimum energy = " << emin << endl;
	cout << "Enter minimum energy = " << emax << endl;
	cout << " "<< endl;
	cout << " "<< endl;	
	

	cout << "AG_add_diff...............................starting"<< endl;		
	if (status == 0)	
		status = AG_add_diff(diffusefilelist, sarfile, edpfile, outfile, emin, emax);
	cout << "AG_add_diff............................... exiting"<< endl;		
	if (status) {
/*		if (status != 105) {
			char * temp = new char[128];
			sprintf(temp, "rm %s", outfile);
			system(temp);
			}*/
		printf("AG_add_diff..................... exiting AG_add_diff ERROR:");		
		fits_report_error(stdout, status);	
		return status;			
		}			
	else {
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_add_diff........... exiting #################\n");
		printf("#################################################################\n\n\n");					
		}			
	
	return status;
}

