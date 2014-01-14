////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       LogFilter pipeline I/O routine
//       GRID event report
//       Release: V0.0 -  16/mar/2005
//       Contributors: A.A., A.G., S.V., S.M.
//       Author: Alessio Trois (IASF-Milano)
//				 Alberto Pellizzoni
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       25/Apr/2005
//                      First release: V1.0
//                      Authors:       Alessio Trois (IASF-Milano)
//                             		   Alberto Pellizzoni(IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

/// #include <alikeLib.h>

#include <cstring>

#include <fstream>
#include <iostream>
#include <string>


/*
#include "fitsio.h"

#include "MathUtils.h"
#include "PilParams.h"
*/
#include "GenmapParams.h"

using namespace std;



static int addfile(fitsfile* iFile, const char* fileList, char* expr, double tmin, double tmax)
{
char buffer[1024];
int status = 0;
FILE *fp = fopen(fileList, "r");
if(!fp) {
	cerr << "Error opening file " << fileList << endl;
	return 104;
	}

int hdutype = 0;
bool noFiles = true;

char name[FLEN_FILENAME];
double t1 = 0, t2 = 0;

while (fgets(buffer , 40960, fp)) {
	sscanf(buffer, "%s %lf %lf", name, &t1, &t2);
	if ( ((t2 > tmin && t2 < tmax)  || (t1 > tmin && t1 < tmax)  || (t1 <= tmin && t2 >= tmax)) && (tmin > obtlimit)) {
		fitsfile *tempFits;
		if (fits_open_file(&tempFits, name, READONLY, &status)) {
			cerr << "Error opening file " << name << endl;
			fclose(fp);
			return status;
			}	
		if (noFiles) {
			fits_copy_file(tempFits, iFile, 1, 1, 1, &status);
			fits_movabs_hdu(iFile, 2, &hdutype, &status);
			fits_select_rows(iFile, iFile, expr, &status);
			fits_close_file(tempFits, &status);
			}
		else {
			fits_movabs_hdu(tempFits, 2, &hdutype, &status);
			fits_select_rows(tempFits, iFile, expr, &status);
			fits_close_file(tempFits, &status);
			}
		noFiles = false;
		}
	}
fclose(fp);
if (noFiles)
	return 1005;
return status;
}





static void EvtExpr(const ThetaGenParams& params, char* expr)
{
/// sprintf(expr, "TIME >= %f && TIME <= %f && ENERGY >= %f && ENERGY <= %f && PH_EARTH > %f && THETA < %f && THETA >= %f", tmin, tmax, emin, emax, albrad, fovradmax, fovradmin);
sprintf(expr, "TIME >= %f && TIME <= %f && ENERGY >= %f && ENERGY <= %f && PH_EARTH > %f && THETA < %f && THETA >= %f", params.tmin, params.tmax, params.emin, params.emax, params.albrad, params.fovradmax, 0.0f);
if ((params.phasecode & 1) == 1) strcat(expr, " && PHASE .NE. 0");
if ((params.phasecode & 2) == 2) strcat(expr, " && PHASE .NE. 1");
if ((params.phasecode & 4) == 4) strcat(expr, " && PHASE .NE. 2");
if ((params.phasecode & 8) == 8) strcat(expr, " && PHASE .NE. 3");
if ((params.phasecode & 16) == 16) strcat(expr, " && PHASE .NE. 4");
if ((params.filtercode & 1) == 1) strcat(expr, " && EVSTATUS .NE. 'L'");
if ((params.filtercode & 2) == 2) strcat(expr, " && EVSTATUS .NE. 'G'");
if ((params.filtercode & 4) == 4) strcat(expr, " && EVSTATUS .NE. 'S'");
}




int countsmalibur(ThetaGenParams& params)
{
	int numcol = 0;
	long nrows = 0; 
	int status = 0;
	double l = 0, b = 0;	
	double x = 0, y = 0;	
	int i = 0, ii = 0;	
	double the = 0;	
	long mxdim = params.mxdim; // dimension (in pixels) of the map
	unsigned short A[mxdim][mxdim];
	double THETA[mxdim][mxdim];
	double THETAVAR[mxdim][mxdim];

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			A[i][ii] = THETA[i][ii] = THETAVAR[i][ii] = 0;
		}
	}	

	double baa = params.ba * DEG2RAD;
	double laa = params.la * DEG2RAD;
	
	int bitpix   =  DOUBLE_IMG; 
	long naxis    =   2;  /* 2-dimensional image                            */    
	long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */		
	
	fitsfile * evtFits;
	char tempname[FLEN_FILENAME];
	strcpy(tempname, tmpnam(NULL));
	if ( fits_create_file(&evtFits, tempname, &status) != 0 ) {
		printf("Errore in apertura file %s\n",tempname);
		return status;
		}	

	char expr[1024];
	EvtExpr(params, expr);
	
	std::cout << std::endl << "AG_thetamapgen....................................adding events files"<< std::endl;

	status = addfile(evtFits, params.evtfile, expr, params.tmin, params.tmax);
	std::cout << "AG_thetamapgen....................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;	
	

	fitsfile * mapFits;
	if ( fits_create_file(&mapFits, params.outfile, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n",params.outfile);
		return status;
		}	

	
	fits_movabs_hdu(evtFits, 2, NULL, &status);	
	fits_get_num_rows(evtFits, &nrows, &status);
	//cout << nrows << endl;
//	double ra[nrows];
//	double dec[nrows];
	double ra, dec, theta;

	string fname(params.outfile);
	fname += ".theta";
	ofstream asciiFile(fname.c_str(), ios::app);
		
	switch (params.projection) {
	case ARC:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			fits_get_colnum(evtFits, 1, "THETA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &theta, NULL, &status);
			//	cout << "RADEC read" << endl;
			/// eulerold(ra, dec, &l, &b, 1);
			Euler(ra, dec, &l, &b, 1);
			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0) {
				the = M_PI;
			} else if (the > 1.0) {
				the = 0.0;
			} else {
				the = acos(the);
			}
			x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
			y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));
	
			i=(int)floor(((-x+(params.mdim/2.))/params.mres));
			ii=(int)floor(((y+(params.mdim/2.))/params.mres));
	
			if (params.inmap(i,ii)) {
					A[ii][i]+=1;
					THETA[ii][i]+=theta;
					//scrivo su un file i theta
					//cout << l << " " << b << " " << laa << " " << baa << endl;
					if(SphDistDeg(l, b, laa, baa) < 2.0 ) {
						asciiFile << theta << " " << ra << " " << dec << endl;
						//cout << theta << " " << ra << " " << dec << endl;
					}
			}
		}
		break;
	case AIT:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			fits_get_colnum(evtFits, 1, "THETA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &theta, NULL, &status);

			/// eulerold(ra, dec, &l, &b, 1);
			Euler(ra, dec, &l, &b, 1);
		    	l *= DEG2RAD;
		    	b *= DEG2RAD;
		    	the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
		    	if (the < -1.0) {
			    the = M_PI;
		    	} else if (the > 1.0) {
			    the = 0.0;
		    	} else {
			    the = acos(the);
		    	}
		    	l=l-laa;

		    	if ( l < M_PI  ) {
		      		l=-l; 
		    	}
		    	else { 
		      		l=2*M_PI -l;
		    	}

		     	x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;         
		     	y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

		    	i=(int)floor(((x+(params.mdim/2.))/params.mres));
		    	ii=(int)floor(((y+(params.mdim/2.))/params.mres));

		   	if (params.inmap(i,ii)) {
			    	A[ii][i]+=1;
				THETA[ii][i]+=theta;
			}
		}
	break;
	}

	asciiFile.close();

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			if(A[i][ii] != 0)
				THETA[i][ii] = THETA[i][ii] / A[i][ii];
			else
				THETA[i][ii] = 0;
		}
	}


	switch (params.projection) {
	    case ARC:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			fits_get_colnum(evtFits, 1, "THETA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &theta, NULL, &status);
			/// eulerold(ra, dec, &l, &b, 1);
			Euler(ra, dec, &l, &b, 1);
			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0) {
				the = M_PI;
			} else if (the > 1.0) {
				the = 0.0;
			} else {
				the = acos(the);
			}
			x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
			y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));
	
			i=(int)floor(((-x+(params.mdim/2.))/params.mres));
			ii=(int)floor(((y+(params.mdim/2.))/params.mres));
	
			if (params.inmap(i,ii)) {
				THETAVAR[ii][i] += pow(theta - THETA[ii][i], 2);
			}
		}
		break;
	    case AIT:
		for (long k = 0; k<nrows; ++k) {
			fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			fits_get_colnum(evtFits, 1, "THETA", &numcol, &status);
			fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &theta, NULL, &status);

			/// eulerold(ra, dec, &l, &b, 1);
			Euler(ra, dec, &l, &b, 1);
		    	l*=DEG2RAD;
		    	b*=DEG2RAD;
		    	the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
		    	if (the < -1.0) {
			    the = M_PI;
		    	} else if (the > 1.0) {
			    the = 0.0;
		    	} else {
			    the = acos(the);
		    	}
		    	l=l-laa;

		    	if ( l < M_PI  ) { 
		      		l=-l; 
		    	}
		    	else { 
		      		l=2*M_PI -l; 
		    	}

		     	x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;         
		     	y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

		    	i=(int)floor(((x+(params.mdim/2.))/params.mres));
		    	ii=(int)floor(((y+(params.mdim/2.))/params.mres));

		   	if (params.inmap(i,ii)) {
				THETAVAR[ii][i] += pow(theta - THETA[ii][i], 2);
			}
		}
		break;
	}

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			if(A[i][ii] != 0)
				THETAVAR[i][ii] = sqrt(THETAVAR[i][ii] / (double)A[i][ii]);
			else
				THETAVAR[i][ii] = 0;
		}
	}	

	/// long nelement =  naxes[0] * naxes[1];
	std::cout<< "creating Theta Map...................................." << std::endl;	
	fits_create_img(mapFits, bitpix, naxis, naxes, &status);
	cout << status << endl;
	std::cout<< "writing Theta Map...................................." << std::endl;
	fits_write_2d_dbl(mapFits, 0, mxdim, mxdim, mxdim, THETA[0], &status);
	cout << status << endl;
	std::cout<< "writing header........................................" << std::endl<< std::endl;	
	params.write_fits_header(mapFits, params.projection, status);
	cout << status << endl;

	
	std::cout<< "creating ThetaVar Map...................................." << std::endl;
 	fits_create_img(mapFits, bitpix, naxis, naxes, &status);
	fits_movabs_hdu(mapFits, 2, 0, &status);
	cout << status << endl;
	std::cout<< "writing ThetaVar Map...................................." << std::endl;
 	fits_write_2d_dbl(mapFits, 0, mxdim, mxdim, mxdim, THETAVAR[0], &status);
	cout << status << endl;
	std::cout<< "writing header........................................" << std::endl<< std::endl;	
 	params.write_fits_header(mapFits, params.projection, status);
	cout << status << endl;
	
	fits_delete_file(evtFits, &status);
// 	cout << status << endl;
	fits_close_file(mapFits, &status);
	cout << status << endl;	
	return status;
	}




int main(int argC, char* argV[])
{
ThetaGenParams params;
if (!params.Load(argC, argV))
	return -1;

int status = 0;


std::cout << " " << std::endl;
std::cout << " " << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "### AG_thetamapgen v.1 - 13/Oct/09 - A.B. (based on ctsmapgen) ##" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << " " << std::endl;


std::cout << "AG_thetamapgen...............................starting"<< std::endl;	
status = countsmalibur(params);
std::cout << "AG_thetamapgen............................... exiting"<< std::endl;


if (status) {
	if (status != 105) {
		remove(params.outfile);
		}
	printf("AG_thetamapgen..................... exiting AG_thetamapgen ERROR:");		
	fits_report_error(stdout, status);	
	return status;
	}			
else {
	printf("\n\n\n###################################################################\n");
	printf("#########  Task AG_thetamapgen........... exiting #################\n");
	printf("#################################################################\n\n\n");					
	}

return status;
}

