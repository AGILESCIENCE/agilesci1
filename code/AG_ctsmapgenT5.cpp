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

#include <sstream>
#include <iostream>
#include <iomanip> // setprecision
#include <fstream>

#include <cstring>

#include <GenmapParams.h>

using namespace std;


string String(const Interval& intv)
{
stringstream str(ios_base::out);
str.precision(6);
str << fixed << intv.Start() << ".." << intv.Stop();
return str.str();
}

string String(const Intervals& intvs)
{
stringstream str(ios_base::out);
const char* sep = "";
for (int i=0; i<intvs.Count(); ++i) {
	str << sep << String(intvs[i]);
	sep = ", ";
	}
return str.str();
}


static string LogEvtString(const Intervals& intvs, double emin, double emax, double albrad, double fovradmax, double fovradmin, int phasecode, int filtercode)
{
if (intvs.Count()<=0)
	return string("");
stringstream str(ios_base::out);
str.precision(6);
if (intvs.Count()==1)
	str << "TIME >= " << fixed << intvs[0].Start() << " && TIME < " << intvs[0].Stop();
else {
	str << "( ";
	const char* sep = "";
	for (int i=0; i<intvs.Count(); ++i) {
		str << sep << " ( TIME >= " << fixed << intvs[i].Start() << " && TIME < " << intvs[i].Stop() << " )";
		sep = " ||";
		}
	str << " )";
	}

str << " && ENERGY >= " << emin;
str << " && ENERGY <= " << emax;
str << " && PH_EARTH > " << albrad;
str << " && THETA < " << fovradmax;
str << " && THETA >= " << fovradmin;

if ((phasecode & 1) == 1)
	str << " && PHASE .NE. 0";
if ((phasecode & 2) == 2)
	str << " && PHASE .NE. 1";
if ((phasecode & 4) == 4)
	str << " && PHASE .NE. 2";
if ((phasecode & 8) == 8)
	str << " && PHASE .NE. 3";
if ((phasecode & 16) == 16)
	str << " && PHASE .NE. 4";

if ((filtercode & 1) == 1)
	str << " && EVSTATUS .NE. 'L'";
if ((filtercode & 2) == 2)
	str << " && EVSTATUS .NE. 'G'";
if ((filtercode & 4) == 4)
	str << " && EVSTATUS .NE. 'S'";

return str.str();
}


static int addfile(fitsfile* iFile, const CtsGenParams& params)
{

std::cout << "params.evtfile: " << params.evtfile << std::endl;

char buffer[1024];
int status = 0;
FILE *fp = fopen(params.evtfile, "r");
if (!fp) {
	cerr << "Error opening file " << params.evtfile << endl;
	return 104;
	}

cout << "Intervals: " << String(params.intervals) << endl;

int hdutype = 0;
bool firstLog = true;
while (fgets(buffer , 40960, fp)) {
	char name[FLEN_FILENAME];
	double t1 = 0, t2 = 0;
	sscanf(buffer, "%s %lf %lf", name, &t1, &t2);
	Interval logIntv(t1, t2);
	Intervals logSel = Intersection(params.intervals, logIntv);

	/**
	if (logSel.Count()) {
		if (logSel.Count()>1)
			cout << name << ": [" << String(logIntv) << "] -> [" << String(logSel) << "]" << endl;
		else if (logIntv==logSel[0])
			cout << name << ": [" << String(logIntv) << "] == [" << String(logSel) << "]" << endl;
		else
			cout << name << ": [" << String(logIntv) << "] -> [" << String(logSel) << "]" << endl;
		}
	*/

	if (logSel.Count())
		cout << "Selecting from: " << name << endl;

	for (int i=0; i<logSel.Count(); ++i) {
		Intervals logSelInt;
		logSelInt.Add(logSel[i]);
		string exprStr = LogEvtString(logSelInt, params.emin, params.emax, params.albrad, params.fovradmax, params.fovradmin, params.phasecode, params.filtercode);
		char expr[1024];
		strcpy(expr, exprStr.c_str());
		fitsfile* tempFits;
		if (firstLog) {
			if (fits_open_file(&tempFits, name, READONLY, &status) != 0 ) {
				cerr << "FITS Error " << status << " opening file " << name << endl;
				return status;
				}	
			fits_copy_file(tempFits, iFile, 1, 1, 1, &status);
			fits_movabs_hdu(iFile, 2, &hdutype, &status);
			fits_select_rows(iFile, iFile, expr, &status);
			}
		else {
			if (fits_open_file(&tempFits, name, READONLY, &status) != 0 ) {
				cerr << "FITS Error " << status << " opening file " << name << endl;
				return status;
				}
			fits_movabs_hdu(tempFits, 2, &hdutype, &status);
			fits_select_rows(tempFits, iFile, expr, &status);
			}
		fits_close_file(tempFits, &status);
		firstLog = false;
		}
	}
fclose(fp);
if (firstLog)
	return 1005;
return status;
}



int countsmalibur(CtsGenParams & params, double timeBin)
{
std::cout << "params.evtfile: " << params.evtfile << std::endl;

/// Creating a temporary fits file to save a first selection of events
int status = 0;
fitsfile* evtFits;
char tempname[FLEN_FILENAME];
strcpy(tempname, tmpnam(NULL));
if (fits_create_file(&evtFits, tempname, &status) != 0 ) {
	cerr << "FITS Error " << status << " opening file " << tempname << endl;
	return status;
	}

const double obtlimit = 104407200.;
if (params.tmin<obtlimit)
	status = 1005;
else
	status = addfile(evtFits, params);
cout << "AG_ctsmapgenT...................................addfile exiting STATUS : "<< status << endl << endl;

ofstream mapText(params.outfile);
streamsize prec = mapText.precision();


int countsCount = int((params.tmax-params.tmin)/timeBin)+1;
int* countsArray = new int[countsCount];
for (int i=0; i<countsCount; ++i)
	countsArray[i] = 0;

fits_movabs_hdu(evtFits, 2, NULL, &status);	
long nrows = 0; 
fits_get_num_rows(evtFits, &nrows, &status);
cout << "Selecting from " << nrows << " rows" << endl;
int selectedRows = 0;
double ra, dec;
for (long k = 0; k<nrows; ++k) {
	int numcol = 0;
	fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
	fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
	fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
	fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
	double l = 0, b = 0;
	Euler(ra, dec, &l, &b, 1);
	if (SphDistDeg(l, b, params.la, params.ba)<=params.mdim/2.0) {
		++selectedRows;
		/// TT_TIME,MJD,L,B,RA,DEC,THETA,PHI,ENERGY,PH_EARTH,EVENT_TYPE
		/// TIME,TIME,L,B,RA,DEC,THETA,PHI,ENERGY,PH_EARTH,EVSTATUS
		double time, theta, phi, energy, phEarth;
		fits_get_colnum(evtFits, 1, "TIME", &numcol, &status);
		fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &time, NULL, &status);
		fits_get_colnum(evtFits, 1, "THETA", &numcol, &status);
		fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &theta, NULL, &status);
		fits_get_colnum(evtFits, 1, "PHI", &numcol, &status);
		fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &phi, NULL, &status);
		fits_get_colnum(evtFits, 1, "ENERGY", &numcol, &status);
		fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &energy, NULL, &status);
		fits_get_colnum(evtFits, 1, "PH_EARTH", &numcol, &status);
		fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &phEarth, NULL, &status);

		fits_get_colnum(evtFits, 1, "EVSTATUS", &numcol, &status);
		char evStatus[128];
		char* evStatusPtr = evStatus;
		fits_read_col(evtFits, TSTRING, numcol, k+1, 1, 1, NULL, &evStatusPtr, NULL, &status);
		mapText.setf(ios::fixed);
		mapText << setprecision(1) << time;
		mapText.unsetf(ios::floatfield);
		mapText.unsetf(ios::fixed);
		mapText <<" "<< setprecision(prec) << l <<" "<< b <<" "<< ra <<" "<< dec <<" "<< theta <<" ";
		mapText << phi <<" "<< energy <<" "<< phEarth <<" "<< evStatus << endl;
		++countsArray[int((time-params.tmin)/timeBin)];
		}
	}
mapText.close();
fits_delete_file(evtFits, &status);

string countsName(params.outfile);
countsName += ".counts";
ofstream countsText(countsName.c_str());
countsText.precision(1);
countsText << fixed;
if (params.tmax==params.tmin+timeBin*(countsCount-1) && countsArray[countsCount-1]==0)
	--countsCount;
for (int i=0; i<countsCount; ++i) {
	Interval timeSlot(params.tmin+timeBin*i, params.tmin+timeBin*(i+1));
	if (Intersection(params.intervals, timeSlot).Count())
		countsText << params.tmin+timeBin*i << " " << params.tmin+timeBin*(i+1) << " " << countsArray[i] << endl;
	}
delete[] countsArray;
countsText.close();

cout << selectedRows << " rows selected" << endl;
return status;
}




int main(int argC, char* argV[])
{


std::cout << " " << std::endl;
std::cout << " " << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "####### AG_ctsmapgenT  V 1.0 - 24/01/12 - A.C., A.T., T.C. ######" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << "#################################################################" << std::endl;
std::cout << " " << std::endl;



CtsGenParams params;
if (!params.Load(argC, argV))
	return -1;

std::cout << "params.evtfile: " << params.evtfile << std::endl;
std::cout << "AG_ctsmapgenT..............................starting"<< std::endl;

double timeBin = 3600; /// zzz ACHTUNG: see the dead code about this

int status = countsmalibur(params, timeBin);
std::cout << "AG_ctsmapgenT.............................. exiting"<< std::endl;

		
if (status) {
	if (status != 105) {
		remove(params.outfile);
		}	
	cerr << "AG_ctsmapgenT.................... exiting AG_ctsmapgen ERROR:";
	fits_report_error(stdout, status);
	return status;
	}			
else {
	printf("\n\n\n#################################################################\n");
	printf("#########  Task AG_ctsmapgenT.......... exiting #################\n");
	printf("#################################################################\n\n\n");					
	}			

return status;
}

