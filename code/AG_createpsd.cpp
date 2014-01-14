////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_createpsd
//       Release: BUILD 0.1 -  29/May/2010
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
//       29/May/2010
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
#include <cstring>
#include <iostream>

#include "pil.h"
#include "fitsio.h"

using namespace std;

	float fitking(float *x, float *params) {
//		TFormula f("kingfunction", "(1. - 1./[1]) * std::pow(1. +  ((x/[0])**2.)/(2.0*[1]), -[1])");
		float f1 = (1. - 1./params[2]) * pow(1.0f +  (x[0]*x[0]/params[1]/params[1])/(2.0f * params[2]), -params[2]);
		float f2 = (1. - 1./params[5]) * pow(1.0f +  (x[0]*x[0]/params[4]/params[4])/(2.0f * params[5]), -params[5]);
		return params[0] * f1 + params[3] * f2;
	}

int AG_createpsd(char *  outfilename, char *  psdtemplate, char *  paramfilename){

	int status;
	fitsfile * psdold;
	fits_open_file(&psdold, psdtemplate, READONLY, &status);
	if (status) {
		cerr << "Cannot open psdtemplate file" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	int bitpix;
	fits_get_img_equivtype(psdold, &bitpix, &status);

// Read rhos and numpsis
	fits_movabs_hdu(psdold, 3, NULL, &status);
	int colnum;

	fits_get_colnum(psdold, CASEINSEN, "RHO", &colnum, &status);
	long nrhos;
	fits_get_coltype(psdold, colnum, NULL, &nrhos, NULL, &status);
	float * rhos = new float[nrhos];
	fits_read_col(psdold, TFLOAT, colnum, 1, 1, nrhos, NULL, rhos, NULL, &status);
	
	fits_get_colnum(psdold, CASEINSEN, "PSI", &colnum, &status);
	long npsis;
	fits_get_coltype(psdold, colnum, NULL, &npsis, NULL, &status);
	
// Open parameter file
	fitsfile * paramfile;
	fits_open_table(&paramfile, paramfilename, READONLY, &status);
	if (status) {
		cerr << "Cannot open parameter file" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	long nrows;
	fits_get_num_rows(paramfile, &nrows, &status);

// Read energies, thetas, and phis from parameter file
	int energycol;
	fits_get_colnum(paramfile, CASEINSEN, "ENERGY", &energycol, &status);
	int thetacol;
	fits_get_colnum(paramfile, CASEINSEN, "THETA", &thetacol, &status);
	int phicol;
	fits_get_colnum(paramfile, CASEINSEN, "PHI", &phicol, &status);
	float energies[1000];
	float thetas[1000];
	float phis[360];
	long nenergies = 0;
	long nthetas = 0;
	long nphis = 0;
	float oldenergy = -1;
	float oldtheta = -1;
	float oldphi = -1;
	for (long row=1; row<nrows; row++) {
		if (oldenergy != -2) {
			fits_read_col(paramfile, TFLOAT, energycol, row, 1, 1, NULL, &energies[nenergies], NULL, &status);
			if (energies[nenergies] < oldenergy) {
				oldenergy = -2;
			} else if (energies[nenergies] > oldenergy) {
				oldenergy = energies[nenergies++];
			}
		}
		if (oldtheta != -2) {
			fits_read_col(paramfile, TFLOAT, thetacol, row, 1, 1, NULL, &thetas[nthetas], NULL, &status);
			if (thetas[nthetas] < oldtheta) {
				oldtheta = -2;
			} else if (thetas[nthetas] > oldtheta) {
				oldtheta = thetas[nthetas++];
			}
		}
		if (oldphi != -2) {
			fits_read_col(paramfile, TFLOAT, phicol, row, 1, 1, NULL, &phis[nphis], NULL, &status);
			if (phis[nphis] < oldphi) {
				oldphi = -2;
			} else if (phis[nphis] > oldphi) {
				oldphi = phis[nphis++];
			}
		}
	}
	long nnewphis = 4 * nphis;
	float  * newphis = new float[nnewphis];
	for (int i=0; i<4; i++) for (int phiind=0; phiind<nphis; phiind++)
		newphis[i*nphis+phiind] = phis[phiind] + 90.0 * i;

// Open output file
	fitsfile * outfile;
	fits_create_file(&outfile, outfilename, &status);
	if (status) {
		cerr << "Cannot create outfile" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	
	long naxes[5];
	naxes[0] = nrhos;
	naxes[1] = npsis;
	naxes[2] = nenergies;
	naxes[3] = nthetas;
	naxes[4] = nnewphis;
	fits_create_img(outfile, bitpix, 5, naxes, &status);
	if (status) {
		cerr << "Cannot create image" << endl;
		fits_report_error(stderr, status);
		return status;
	}

// Fill psd matrix
	int normcol;
	fits_get_colnum(paramfile, CASEINSEN, "NORM", &normcol, &status);
	int gammacol;
	fits_get_colnum(paramfile, CASEINSEN, "GAMMA", &gammacol, &status);
	int angcol;
	fits_get_colnum(paramfile, CASEINSEN, "ANG", &angcol, &status);
	if (status) {
		cerr << "Cannot find parameter table column numbers" << endl;
		fits_report_error(stderr, status);
		return status;
	}

	float params[6]={0.0,0.0,0.0,0.0,1.0,1.0};
	long row=1;
	long fpixel[5];
	for (fpixel[2]=1; fpixel[2]<=nenergies; fpixel[2]++)
		for (fpixel[3]=1; fpixel[3]<=nthetas; fpixel[3]++)
			for (int phiind=1; phiind<=nphis; phiind++){
				fits_read_col(paramfile, TFLOAT, normcol, row, 1, 1, NULL, &params[0], NULL, &status);
				fits_read_col(paramfile, TFLOAT, angcol, row, 1, 1, NULL, &params[1], NULL, &status);
				fits_read_col(paramfile, TFLOAT, gammacol, row++, 1, 1, NULL, &params[2], NULL, &status);
				for (int rhoind=0; rhoind<nrhos; rhoind++){
					float king = fitking(&rhos[rhoind],params);
					fpixel[0] = rhoind+1;
					for (fpixel[1] = 1; fpixel[1]<=npsis; fpixel[1]++)
						for(int quad=0; quad<4; quad++) {
							fpixel[4] = quad*nphis+phiind;
							fits_write_pix(outfile, TFLOAT, fpixel, 1, &king, &status);
					}
				}
	}
	if (status) {
		cerr << "Problems filling PSD matrix" << endl;
		fits_report_error(stderr, status);
		return status;
	}

	delete [] rhos;
	fits_close_file(paramfile, &status);
	if (status) {
		cerr << "Problems closing parameter file" << endl;
		fits_report_error(stderr, status);
		return status;
	}

// Copy primary image keywords
	fits_movabs_hdu(psdold, 1, NULL, &status);
	if (status) {
		cerr << "Problems moving to primary HDU in PSD template file" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	int numkeys;
	fits_get_hdrspace(psdold, &numkeys, NULL, &status);
	if (status) {
		cerr << "Problems reading number of keywords in primary HDU" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	for (int key=1; key<=numkeys; key++){
		char keyname[FLEN_FILENAME];
		char value[FLEN_FILENAME];
		char comment[FLEN_FILENAME];
		char card[FLEN_FILENAME];
		fits_read_keyn(psdold, key, keyname, value, comment, &status);
		if (strncmp(keyname, "NAXIS", 5)) {
			fits_read_record(psdold, key, card, &status);
			fits_update_card(outfile, keyname, card, &status);
		}
	}
	fits_write_chksum(outfile, &status);
	if (status) {
		cerr << "Cannot copy primary header keywords" << endl;
		fits_report_error(stderr, status);
		return status;
	}

// Copy template file HDUs to output file

	fits_copy_file(psdold, outfile, 0, 0, 1, &status);
	fits_close_file(psdold, &status);
	if (status) {
		cerr << "Cannot copy extensions" << endl;
		fits_report_error(stderr, status);
		return status;
	}

// Overwrite energies, thetas, and phis in output file

	int typecode;

	fits_get_colnum(outfile, CASEINSEN, "ENERGY", &colnum, &status);
	fits_modify_vector_len(outfile, colnum, nenergies, &status);
	fits_write_col(outfile, TFLOAT, colnum, 1, 1, nenergies, energies, &status);
	
	fits_get_colnum(outfile, CASEINSEN, "POLAR_ANGLE", &colnum, &status);
	fits_modify_vector_len(outfile, colnum, nthetas, &status);
	fits_write_col(outfile, TFLOAT, colnum, 1, 1, nthetas, thetas, &status);

	fits_get_colnum(outfile, CASEINSEN, "AZIMUTH_ANGLE", &colnum, &status);
	fits_modify_vector_len(outfile, colnum, nnewphis, &status);
	fits_write_col(outfile, TFLOAT, colnum, 1, 1, nnewphis, newphis, &status);

	if (status) {
		cerr << "Cannot overwrite axes" << endl;
		fits_report_error(stderr, status);
		return status;
	}
	fits_write_chksum(outfile, &status);

	delete [] newphis;
	fits_close_file(outfile, &status);
	fits_report_error(stderr, status);
	return status;
}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char * outfilename = new char[FLEN_FILENAME];
	char * psdtemplate  = new char[FLEN_FILENAME];
	char * paramfilename  = new char[FLEN_FILENAME];
	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("outfile", outfilename);
	status = PILGetString("psdtemplate", psdtemplate);
	status = PILGetString("paramfile", paramfilename);

	status = PILClose(status);

	if (!status) {
	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_createpsd.cpp v.0.1 -29/05/10 - A.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Output file name : " <<  outfilename << endl;
	cout << "PSD template file name : "<< psdtemplate << endl;
	cout << "Parameter file name : "<< paramfilename << endl;
	

	cout << "AG_createpsd...............................starting"<< endl;		
	status = AG_createpsd(outfilename, psdtemplate, paramfilename);
	cout << "AG_createpsd............................... exiting"<< endl;		
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_createpsd........... exiting #################\n");
		printf("#################################################################\n\n\n");					
	
	} else
		cerr << "PIL error" << endl;
	delete[] outfilename;
	delete[] psdtemplate;
	delete[] paramfilename;

	return status;
}

