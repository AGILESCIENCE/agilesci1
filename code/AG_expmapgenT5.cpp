////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       excalibur
//       Release: V2.5 - 27/Aug/2011
//		Contributors: 
//       Author: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
//               Andrea Bulgarelli (IASF-Bologna), Tomaso Contessi (Nuove Idee sas)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//		 V2.5 -  27/Aug/2011
//		 V2.2 -  31/Jul/2009
//		 V2.1 -  23/Jan/2009
//       V1.0 -  8/Dec/2005
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "pil.h"

#include "wcshdr.h"
#include "wcsmath.h"
#include "wcstrig.h"
#include "sph.h"

#include "CalibUtils.h"
#include "FitsUtils.h"
#include "GenmapParams.h"



using namespace std;



/// s_localFileName is used first to copy the remote log files, and later to store the selection
char s_localFileName[FLEN_FILENAME];
char s_selectionFileName[FLEN_FILENAME];	/// Destination file for the selected events



static string String(const Interval& intv)
{
stringstream str(ios_base::out);
str.precision(6);
str << fixed << intv.Start() << ".." << intv.Stop();
return str.str();
}

static string String(const Intervals& intvs)
{
stringstream str(ios_base::out);
const char* sep = "";
for (int i=0; i<intvs.Count(); ++i) {
	str << sep << String(intvs[i]);
	sep = ", ";
	}
return str.str();
}



static string LogExprTimeString(const Intervals& intvs)
{
if (intvs.Count()<=0)
	return string("");
stringstream str(ios_base::out);
str.precision(6);
if (intvs.Count()==1)
	str << "TIME >= " << fixed << intvs[0].Start() << " && TIME < " << intvs[0].Stop();
else {
	str << "(";
	const char* sep = "";
	for (int i=0; i<intvs.Count(); ++i) {
		str << sep << "(TIME >= " << fixed << intvs[i].Start() << " && TIME < " << intvs[i].Stop() << ")";
		sep = " || ";
		}
	str << ")";
	}
return str.str();
}


static string LogExprString(const Intervals& intvs, int phasecode, int timeStep)
{
if (intvs.Count()<=0)
	return string("");
stringstream str(ios_base::out);

str << LogExprTimeString(intvs);

str << " && LIVETIME > 0 && LOG_STATUS == 0 && MODE == 2";

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
str << " && ((#ROW == 1) || (#ROW == (#ROW/" << timeStep << ") *" << timeStep << "))";
return str.str();
}


const double c_binFactor = 0.0003046174197867085688996857673060958405;
const double c_angleFactor = 0.0174532925199432954743716805978692718782;

double Area(double xbin, double ybin, double theta, int projection)
{
if (projection==ARC)
	return c_binFactor * xbin * ybin * Sinaa(c_angleFactor * theta);
else /// AIT
	return c_binFactor * xbin * ybin;
}

/**
static void ReadFitsCol(FitsFile& file, double* array, const char* colName, long rowOffs, long rowCount)
{
int colNum;
fits_get_colnum(file, CASESEN, const_cast<char*>(colName), &colNum, file);
fits_read_col(file, TDOUBLE, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, file);
}

static void ReadFitsCol(fitsfile* file, short* array, const char* colName, long rowOffs, long rowCount)
{
int colNum;
fits_get_colnum(file, CASESEN, const_cast<char*>(colName), &colNum, file);
fits_read_col(file, TSHORT, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, file);
}
*/

static void ReadFitsCol(fitsfile* file, double* array, const char* colName, long rowOffs, long rowCount, int* status)
{
int colNum;
fits_get_colnum(file, CASESEN, const_cast<char*>(colName), &colNum, status);
fits_read_col(file, TDOUBLE, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}

static void ReadFitsCol(fitsfile* file, short* array, const char* colName, long rowOffs, long rowCount, int* status)
{
int colNum;
fits_get_colnum(file, CASESEN, const_cast<char*>(colName), &colNum, status);
fits_read_col(file, TSHORT, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}



/// Some variable were accessed but not used in previous versions.
/// Their definition and usage has been put under ifdef __UNUSED__VARS__

/// #define __UNUSED__VARS__

#ifdef __UNUSED__VARS__
#define ALLOCATE_ARRAYS(count) \
	double* q1 = new double[count];\
	double* q2 = new double[count];\
	double* q3 = new double[count];\
	double* q4 = new double[count];\
	short*  mode = new short[count];\
	short*  log_status = new short[count]; \
	double* ra_y = new double[count];\
	double* dec_y = new double[count];\
	double* psi = new double[count];\
	double* gp = new double[count];\
	double* livetime = new double[count];\
	double* earth_ra = new double[count];\
	double* earth_dec = new double[count];\
	short*  phase = new short[count];\
	long* change = new long[count]
#else
#define ALLOCATE_ARRAYS(count) \
	double* ra_y = new double[count];\
	double* dec_y = new double[count];\
	double* psi = new double[count];\
	double* gp = new double[count];\
	double* livetime = new double[count];\
	double* earth_ra = new double[count];\
	double* earth_dec = new double[count];\
	short*  phase = new short[count];\
	long* change = new long[count]
#endif

#ifdef __UNUSED__VARS__
#define DELETE_ARRAYS() \
	delete[] q1;\
	delete[] q2;\
	delete[] q3;\
	delete[] q4;\
	delete[] mode;\
	delete[] log_status; \
	delete[] ra_y;\
	delete[] dec_y;\
	delete[] psi;\
	delete[] gp;\
	delete[] livetime;\
	delete[] earth_ra;\
	delete[] earth_dec;\
	delete[] phase;\
	delete[] change
#else
#define DELETE_ARRAYS() \
	delete[] ra_y;\
	delete[] dec_y;\
	delete[] psi;\
	delete[] gp;\
	delete[] livetime;\
	delete[] earth_ra;\
	delete[] earth_dec;\
	delete[] phase;\
	delete[] change
#endif


static int pil_chridx(char *p, char c)
{
for (int l=0; p[l]; l++)
	if (c == p[l])
		return l;
return -1;
}


static int pil_curly_append(char **p, char *s, int nchars)
{
	int	l, nl;
	char	*np;
	
	if (NULL == p)  return(PIL_NUL_PTR);
	if (nchars < 0) return(PIL_BAD_ARG);
	if (0 == nchars)
	{ if (NULL != *p) return(PIL_OK); }
	else
	{ if (NULL == s)  return(PIL_NUL_PTR); }
	
	if (NULL == *p)  { l = 0; }
	else  { l = strlen(*p); }
	
	nl = nchars + l + 1;						/* + 1 - for EOS */
	
	if (NULL == *p)  { np = (char *)PIL_malloc(nl); }
	else  { np = (char *)PIL_realloc((void *)*p, nl); }
	if (NULL == np)
	{ if (NULL != *p)  { PIL_free(*p); *p = NULL; }		/* on error deallocate string */
		return(PIL_NO_MEM);
	}
	
	if (nchars > 0) memcpy(np + l, s, nchars);			/* append new string */
	np[nl - 1] = 0;						/* signal EOS */
	
	*p = np;
	return(PIL_OK);
}


static int pil_curly_expand(char *src, char **dst)
{ char		ev[PIL_CURLYSIZE];
	int		sl, l, l2, l3, tl, r;
	char		*p;
	
	
	if (NULL == src) return(PIL_NUL_PTR);
	if (NULL == dst) return(PIL_NUL_PTR);
	*dst = NULL;
	
	sl = tl = 0;
	for (;;)
    { if (-1 == (l = pil_chridx(src + sl + tl, '$'))) break;	/* no '$' from current pos, so terminate */
		if ('{' != src[sl + tl + l + 1])
        { tl += l + 1;
			continue;						/* '$' not followed by '{' */
        }
		
		l += tl;							/* total number of bytes _BEFORE_ '${...}' */
		if (-1 == (l2 = pil_chridx(src + sl + l + 2, '}'))) break;  /* no closing '}' */
		if (PIL_OK != (r = pil_curly_append(dst, src + sl, l))) return(r);  /* copy part _BEFORE_ '$' */
		sl += l + 2;						/* position on 1st char after '${' */
		
		if (l2 > 0)						/* if non-empty env.var name */
        { memcpy(ev, src + sl, l2);				/* then copy its name to the temporary buffer */
			ev[l2] = 0;						/* signal end of string */
			p = getenv(ev);					/* get environment variable */
			l3 = 0;
			if (NULL != p)  l3 = strlen(p);			/* ... and compute length of that variable (0 - if notdef) */
			if (PIL_OK != (r = pil_curly_append(dst, p, l3))) return(r);  /* copy value of environment variable */
		}
		sl += l2 + 1;						/* env.var.name + '}' */
		tl = 0;
    }
	
	return(pil_curly_append(dst, src + sl, strlen(src + sl)));	/* copy remaining part */
}



static bool CopyFile(const char* iFileName, const char* oFileName)
{
const int MAX_BUFF_SIZE = 65536;
static char s_fileBuff[MAX_BUFF_SIZE];

ifstream iFile(iFileName, ios::in|ios::binary);
if (!iFile.is_open()) {
	cerr << "CopyFile ERROR: cannot open " << iFileName << " for reading" << endl;
	return false;
	}
ofstream oFile(oFileName, ios::out|ios::binary);
if (!oFile.is_open()) {
	cerr << "CopyFile ERROR: cannot open " << oFileName << " for writing" << endl;
	return false;
	}

iFile.seekg(0, ios::end);
long fileSize = iFile.tellg();
iFile.seekg(0, ios::beg);

cout << "Copying " << fileSize << " bytes from " << iFileName << " to " << oFileName << "." << flush;
int progress = 0;
while (fileSize>0) {
	int transSize = fileSize>MAX_BUFF_SIZE?MAX_BUFF_SIZE:fileSize;
	iFile.read(s_fileBuff, transSize);
	oFile.write(s_fileBuff, transSize);
	fileSize -= transSize;
	++progress;
	if (progress%32==0)
		cout << "." << flush;
	}
cout << "Done" << endl;
return true;
}


static int MakeLogSelection(const char* fileList, const Intervals& selection, int phasecode, int timeStep)
{
FILE *fp = fopen(fileList, "r");
if (!fp) {
	cerr << "Error opening file " << fileList << endl;
	return -100; /// already printed error
	}

cout << "MakeLogSelection: Intervals: " << String(selection) << endl;

char scratchFileName[FLEN_FILENAME];
tmpnam(scratchFileName);
strcat(scratchFileName, ".gz");

int skippedFiles = 0;

int status = 0;
/// fitsfile* tempFits = 0;
FitsFile tempFits;


char buffer[1024];
while (fgets(buffer, 1023, fp) && !status) {
	char name[FLEN_FILENAME];
	double t1=0, t2=0;
	sscanf(buffer, "%s %lf %lf", name, &t1, &t2);

	Interval logIntv(t1, t2);
	Intervals logSel = Intersection(selection, logIntv);

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

	char* nname = 0;
	pil_curly_expand(name, &nname);
	if (logSel.Count()) {
		if (skippedFiles) {
			cout << "Skipped " << skippedFiles << " file: " << endl;
			skippedFiles = 0;
			}
		cout << "Selecting from: " << nname << endl;
		}
	else
		++skippedFiles;
		/// cout << "Skipping file: " << nname << endl;

	for (int i=0; i<logSel.Count() && !status; ++i) {

		
		Intervals logSelInt;
		logSelInt.Add(logSel[i]);
		string exprStr = LogExprString(logSelInt, phasecode, timeStep);
		char expr[1024];
		strcpy(expr, exprStr.c_str());

		if (!CopyFile(nname, scratchFileName)) {
			PIL_free(nname);
			fclose(fp);
			cerr << "Error copying " << nname << " to " << scratchFileName << endl;
			return -100; /// already printed error
			}

/**
		cerr << "Unzipping " << copytempnamegz << " to " << copytempname << endl;
		cerr << gunzip << endl;
		system(gunzip);
		cerr << "Done" << endl;
*/

		/// fitsfile *iFile;
		/// if (fits_open_file(&iFile, scratchFileName, READONLY, &status)) {
		FitsFile iFile;
		if (!iFile.Open(scratchFileName)) {
			cerr << "Failed opening: " << scratchFileName << ", copy of " << nname << endl;
			PIL_free(nname);
			fclose(fp);
			return -100; /// already printed error
			}
		int hdutype = 0;
		if (tempFits.File()) {
			iFile.MoveAbsHDU(2); /// fits_movabs_hdu(iFile, 2, &hdutype, &status);
			fits_select_rows(iFile, tempFits, expr, &status);
			}
		else {
			/// if (fits_create_file(&tempFits, s_selectionFileName, &status)) {
			if (!tempFits.Create(s_selectionFileName)) {
				cerr << "Error creating the events selection file " << s_selectionFileName << endl;
				PIL_free(nname);
				fclose(fp);
				return -100; /// already printed error
				}

			fits_copy_file(iFile, tempFits, 1, 1, 1, &status);
			tempFits.MoveAbsHDU(2); /// fits_movabs_hdu(tempFits, 2, &hdutype, &status);
			fits_select_rows(tempFits, tempFits, expr, &status);

			if (status==0) {

/**
				fitsfile* templateFile;
				if (fits_create_file(&templateFile, s_localFileName, &status))
					cerr << "Error creating the template local file " << s_localFileName << endl;
				else {
					fits_copy_file(tempFits, templateFile, 1, 1, 1, &status);
					char timeSel[12] = "TIME < 0";
					fits_select_rows(templateFile, templateFile, timeSel, &status); /// zzz Just to delete all rows (find a better way)
					fits_close_file(templateFile, &status);
					}
*/
				FitsFile templateFile;
				if (!templateFile.Create(s_localFileName))
					cerr << "Error creating the template local file " << s_localFileName << endl;
				else {
					fits_copy_file(tempFits, templateFile, 1, 1, 1, &status);
					char timeSel[12] = "TIME < 0";
					fits_select_rows(templateFile, templateFile, timeSel, templateFile); /// zzz Just to delete all rows (find a better way)
					}


				}
			}
		/// fits_close_file(iFile, &status);
		}
	PIL_free(nname);
	}
fclose(fp);
if (skippedFiles)
	cout << "Skipped " << skippedFiles << " file: " << endl;

if (tempFits.File()) {
	/// fits_close_file(tempFits, &status);
	remove(scratchFileName);
	}
else if (!status)
	status = -118;	/// No events selected

return status;
}




bool EvalExposure(fitsfile* localFits, const AeffGridAverage& raeff, ExpGenParamsT& params, const Intervals& intvs, double* resultingExp)
{
double lp = 0, bp = 0;
double learth, bearth;
double lp0 = 0, bp0 = 0;
double x = 0, y = 0;
double theta = 0, phi = 0, phi2 = 0;
double eul[5];


double theta2;
double lng, lat;
double area = 0;
int ait;

struct prjprm prj;
prjini(&prj);

double mdim = params["mdim"];

switch (params.projection) {
case ARC:
	x = -mdim;
	y = -mdim;
	ait = 0;
	theta2 = 90.0-sqrt(x*x+y*y);
	phi2 = atan2d(-y, -x);
	eul[0] = params["la"];
	eul[1] = 90.0-double(params["ba"]);
	eul[2] = params["lonpole"];
	eul[3] = cosd(eul[1]);
	eul[4] = sind(eul[1]);	
		
	sphx2s(eul, 1, 1, 0, 0, &phi2, &theta2, &lng, &lat);
	area = Area(mdim, mdim, 90-theta2, params.projection);
	break;
case AIT:
	prj.r0 = RAD2DEG;
	prj.phi0 = params["la"];
	prj.theta0 = params["ba"];
	aitset(&prj);
	x = mdim;
	y = -mdim;
	aitx2s(&prj, 1, 1, 0, 0, &x, &y, &lng, &lat, &ait);
	area = Area(mdim, mdim, 0, params.projection);
	}

long n = 0;
double time = 0;

/// Copying a selection of the selection to the local file

int status = 0;
fitsfile* tempFits;
if (fits_open_file(&tempFits, s_localFileName, READWRITE, &status)) {
	cerr << "Error opening the events template file " << s_localFileName << endl;
	return status;
	}
int hdutype = 0;
fits_movabs_hdu(tempFits, 2, &hdutype, &status);
char timeSel[12] = "TIME < 0";
fits_select_rows(tempFits, tempFits, timeSel, &status); /// zzz Just to delete all rows (find a better way)

string exprStr = LogExprTimeString(intvs);
char expr[1024];
strcpy(expr, exprStr.c_str());

fits_select_rows(localFits, tempFits, expr, &status);
if (status) {
	cerr << "Error adding rows to the second selection file " << s_localFileName << endl;
	return 105;
	}
long allnrows;
fits_get_num_rows(tempFits, &allnrows, &status);

long rowblockincrement = 65536;
fits_get_rowsize(tempFits, &rowblockincrement, &status);
cout << "rowblockincrement=" << rowblockincrement << endl;

ALLOCATE_ARRAYS(rowblockincrement);

double A = 0.0; /// The resulting value
/**
double y_tol = params["y_tol"];
double earth_tol = params["earth_tol"];
double fovrad = params["fovrad"];
double fovradmin = params["fovradmin"];
*/

double timestep = params["timestep"];

for (long rowblockzero=0; rowblockzero < allnrows; rowblockzero += rowblockincrement) {
	long nrows = (rowblockzero + rowblockincrement < allnrows) ? rowblockincrement : (allnrows - rowblockzero);
#ifdef __UNUSED__VARS__
	ReadFitsCol(tempFits, q1, "Q1", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, q2, "Q2", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, q3, "Q3", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, q4, "Q4", rowblockzero, nrows, &status);
#endif
	ReadFitsCol(tempFits, ra_y, "ATTITUDE_RA_Y", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, dec_y, "ATTITUDE_DEC_Y", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, livetime, "LIVETIME", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, phase, "PHASE", rowblockzero, nrows, &status);
#ifdef __UNUSED__VARS__
	ReadFitsCol(tempFits, mode, "MODE", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, log_status, "LOG_STATUS", rowblockzero, nrows, &status);
#endif
	ReadFitsCol(tempFits, earth_ra, "EARTH_RA", rowblockzero, nrows, &status);
	ReadFitsCol(tempFits, earth_dec, "EARTH_DEC", rowblockzero, nrows, &status);
	
	long count = 0;
	double earth_ra0 = earth_ra[0], earth_dec0 = earth_dec[0];
	double ra_y0 = ra_y[0], dec_y0 = dec_y[0];
	for (long k = 1; k<nrows; ++k) {
		if ((phase[k-1] != phase[k])
				|| params.YTolTest(ra_y[k-1], dec_y[k-1], ra_y0, dec_y0)
				|| params.EarthTolTest(earth_ra[k-1], earth_dec[k-1], earth_ra0, earth_dec0)) {
			change[count++] = k;
			earth_ra0 = earth_ra[k-1];
			earth_dec0 = earth_dec[k-1];
			ra_y0 = ra_y[k-1];
			dec_y0 = dec_y[k-1];
			}
		}
	if (count == 0)
		change[0] = nrows;

	long lowrow = 0;
	long highrow = 0;

	for (long k = 0; k<=count; ++k) {
		time = 0;
		if ( k == 0 ) {
			lowrow = 0;
			highrow = change[0];
			}
		else if ( k == count) {
			lowrow = change[k-1];
			highrow = nrows;
			}
		else {
			lowrow = change[k-1];
			highrow = change[k];
			}
	
		for (n = lowrow; n<highrow; n++)
			time += livetime[n] * timestep;

		Euler(ra_y[lowrow], dec_y[lowrow], &lp, &bp, 1);
		Euler(earth_ra[lowrow], earth_dec[lowrow], &learth, &bearth, 1);
		if (k==0 && rowblockzero==0) {
			lp0 = lp;
			bp0 = bp;
			}

		if (ait == 0) {
			theta = SphDistDeg(lng, lat, lp, bp);
			phi = 0;
			if (params.FovTest(theta) && params.AlbTest(lng, lat, learth, bearth))
				A += 1e-3*time*(raeff.AvgVal(theta, phi))*area;
			}
		}
	}

*resultingExp = A;

fits_close_file(tempFits, &status);
DELETE_ARRAYS();
return status==0;
}


class AppScreen
{
public:
	AppScreen()
	{
	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "### AG_expmapgenT V1.0 - 25/01/2012 - A.C., A.T., A.B., T.C. ####"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << endl;
	}

	~AppScreen()
	{
	cout << endl << endl << endl;
	cout << "#################################################################"<< endl;
	cout << "##########  Task AG_expmapgenT........... exiting ###############"<< endl;
	cout << "#################################################################"<< endl << endl;
	}
};



static void PrintVector(const VecF& arr)
{
cout << "[" << arr.Size() << "] ";
for (int i=0; i<arr.Size(); ++i)
	cout << " " << arr[i];
cout << endl;
}


static void TestEdp(const char* fileName)
{
EdpGrid edpGrid;
int status = edpGrid.Read(fileName);

const VecF& trueEnes = edpGrid.TrueEnergies();
const VecF& obsEnes = edpGrid.ObsEnergies();
const VecF& thetas = edpGrid.Thetas();
const VecF& phis = edpGrid.Phis();
cout << "TrueEnergies: "; PrintVector(trueEnes);
cout << "ObsEnergies: "; PrintVector(obsEnes);
cout << "Thetas: "; PrintVector(thetas);
cout << "Phis: "; PrintVector(phis);


EdpGridClass edpGridOld;
status = edpGridOld.Read(fileName);

const Mat4F& edpGridVals = edpGrid.Values();
cout << "edpGridVals.Dim(0): " << edpGridVals.Dim(0) << endl;
cout << "edpGridVals.Dim(1): " << edpGridVals.Dim(1) << endl;
cout << "edpGridVals.Dim(2): " << edpGridVals.Dim(2) << endl;
cout << "edpGridVals.Dim(3): " << edpGridVals.Dim(3) << endl;

cout << "edpGridOld.DimOf(0): " << edpGridOld.DimOf(0) << endl;
cout << "edpGridOld.DimOf(1): " << edpGridOld.DimOf(1) << endl;
cout << "edpGridOld.DimOf(2): " << edpGridOld.DimOf(2) << endl;
cout << "edpGridOld.DimOf(3): " << edpGridOld.DimOf(3) << endl;


/// cout << "edpGrid.Val(12,2,2,317): " << edpGrid.Val(12,2,2,317) << endl;
/// cout << "edpGridOld.Val(12,2,2,317): " << edpGridOld.Val(12,2,2,317) << endl;


float D = -2;
int diffCountMat = 0;
int diffCountVal = 0;

for (int i0 = 0; i0<trueEnes.Size(); ++i0)
	for (int i1 = 0; i1<obsEnes.Size(); ++i1)
		for (int i2 = 0; i2<thetas.Size(); ++i2)
			for (int i3 = 0; i3<phis.Size(); ++i3) {
/// cerr << i0 << " " << i1 << " " << i2 << " " << i3 << endl;
/// cerr << trueEnes[i0]+D << " " << obsEnes[i1]+D << " " << thetas[i2]+D << " " << phis[i3]+D << endl;

					if (edpGridOld.ValOf(i3,i2,i1,i0) != edpGridVals(i3,i2,i1,i0))
						++diffCountMat;
					if (edpGridOld.Val(trueEnes[i0]+D,obsEnes[i1]+D,thetas[i2]+D,phis[i3]+D) != edpGrid.Val(trueEnes[i0]+D,obsEnes[i1]+D,thetas[i2]+D,phis[i3]+D))
						++diffCountVal;
				}
cout << "diffCountMat = " << diffCountMat << endl;
cout << "diffCountVal = " << diffCountVal << endl;

}


static void TestPsf(const char* fileName)
{
PsfGrid psfGrid;
int status = psfGrid.Read(fileName);

const VecF& rhos = psfGrid.Rhos();
const VecF& psis = psfGrid.Psis();
const VecF& thetas = psfGrid.Thetas();
const VecF& phis = psfGrid.Phis();
const VecF& energies = psfGrid.Energies();

cout << "Rhos: "; PrintVector(rhos);
cout << "Psis: "; PrintVector(psis);
cout << "Energies: "; PrintVector(energies);
cout << "Thetas: "; PrintVector(thetas);
cout << "Phis: "; PrintVector(phis);


PsfGridClass psfGridOld;
status = psfGridOld.Read(fileName);

const Mat5F& psfGridVals = psfGrid.Values();
cout << "psfGridVals.Dim(0): " << psfGridVals.Dim(0) << endl;
cout << "psfGridVals.Dim(1): " << psfGridVals.Dim(1) << endl;
cout << "psfGridVals.Dim(2): " << psfGridVals.Dim(2) << endl;
cout << "psfGridVals.Dim(3): " << psfGridVals.Dim(3) << endl;
cout << "psfGridVals.Dim(4): " << psfGridVals.Dim(4) << endl;

cout << "psfGridOld.DimOf(0): " << psfGridOld.DimOf(0) << endl;
cout << "psfGridOld.DimOf(1): " << psfGridOld.DimOf(1) << endl;
cout << "psfGridOld.DimOf(2): " << psfGridOld.DimOf(2) << endl;
cout << "psfGridOld.DimOf(3): " << psfGridOld.DimOf(3) << endl;
cout << "psfGridOld.DimOf(4): " << psfGridOld.DimOf(4) << endl;

/// cout << "psfGrid.Val(12,2,2,317): " << psfGrid.Val(12,2,2,317) << endl;
/// cout << "psfGridOld.Val(12,2,2,317): " << psfGridOld.Val(12,2,2,317) << endl;

int diffCountMat = 0;
int diffCountVal = 0;
for (int i0 = 0; i0<rhos.Size(); ++i0)
	for (int i1 = 0; i1<psis.Size(); ++i1)
		for (int i2 = 0; i2<energies.Size(); ++i2)
			for (int i3 = 0; i3<thetas.Size(); ++i3)
				for (int i4 = 0; i4<phis.Size(); ++i4) {

/// cerr << i0 << " " << i1 << " " << i2 << " " << i3 << i4 << endl;
/// cerr << rhos[i0] << " " << psis[i1] << " " << thetas[i2] << " " << phis[i3]  << " " << energies[i4] << endl;

					if (psfGridOld.ValOf(i4,i3,i2,i1,i0) != psfGridVals(i4,i3,i2,i1,i0))
						++diffCountMat;
///					if (psfGridOld.Val(rhos[i0]+D,psis[i1]+D,energies[i2]+D,thetas[i3]+D,phis[i4]+D) != psfGrid.Val(rhos[i0]+D,psis[i1]+D,energies[i2]+D,thetas[i3]+D,phis[i4]+D))
///						++diffCountVal;
					}
cout << "diffCountMat = " << diffCountMat << endl;
cout << "diffCountVal = " << diffCountVal << endl;
}



int main(int argc,char **argv)
{
ExpGenParamsT params;
if (!params.Load(argc, argv))
	return -1;


// #define JUST_TESTING

#ifdef JUST_TESTING

/// TestEdp(params["raeffFileName"]);

TestPsf(params["raeffFileName"]);

return 0;

#endif

Intervals intvs;
double tmin = params["tmin"];
double tmax = params["tmax"];
const char* intFileName = params["timelist"];

if (strcmp(intFileName, "None")) {
	intvs = ReadIntervals(intFileName);
	tmin = intvs.Min();
	tmax = intvs.Max();
	params["tmin"] = tmin;
	params["tmax"] = tmax;
	}
else {
	Interval intv(tmin, tmax);
	intvs.Add(intv);
	}

params.Print();
if (intvs.Count()>1) {
	cout << intvs.Count() << " intervals:" << endl;
	for (int i=0; i<intvs.Count(); ++i)
		cout << "   " << String(intvs[i]) << endl;
	}

cout << "AG_expmapgenT......................selecting the events"<< endl;

tmpnam(s_selectionFileName);
strcat(s_selectionFileName, ".gz");
tmpnam(s_localFileName);

const char* outfile = params["outfile"];

const char* logfile = params["logfile"];
if (logfile && logfile[0]=='@')
	++logfile;

int status = MakeLogSelection(logfile, intvs, params["phasecode"], params["timestep"]);

if (status==-118)
	cout << endl << "AG_expmapgenT......................no matching events found"<< endl;
else if (status==-100)
	cout << endl << "AG_expmapgenT......................selection failed"<< endl;
else if (status!=0)
	cout << endl << "AG_expmapgenT......................selection failed"<< endl;
else {

	double deltaT = params["timeslot"];

	cout << "AG_expmapgenT......................evaluating the exposure"<< endl;

	double emin = params["emin"];
	double emax = params["emax"];
	double index = params["spectral_index"];
	AeffGridAverage raeff(params["raeffFileName"], emin, emax, index);

	fitsfile* localFits;
	if (fits_open_file(&localFits, s_selectionFileName, READONLY, &status)) {
		cerr << "ERROR opening local selection file " << s_selectionFileName << endl;
		return -1;
		}
	int hdutype = 0;
	fits_movabs_hdu(localFits, 2, &hdutype, &status);


	ofstream expText(outfile);
	expText.setf(ios::fixed);
	double beginTime = tmin;
	double endTime = beginTime+deltaT;
	if (endTime>tmax)
		endTime = tmax;
	cout << "***** " << beginTime << " " << endTime << " " << deltaT << endl << endl;
	double totalExposure = 0;
	Interval timeSlot;
	do {
		timeSlot.Set(beginTime, endTime);
		Intervals intervalSlots = Intersection(intvs, timeSlot);
		if (intervalSlots.Count()) {

			cout << "Selected intervals" << endl;
			for (int i=0; i<intervalSlots.Count(); ++i)
				cout << "   " << intervalSlots[i].Start() << " " << intervalSlots[i].Stop() << endl;
			double exp;
			if (EvalExposure(localFits, raeff, params, intervalSlots, &exp)) {
				expText << setprecision(1);
				expText << beginTime << " " << endTime << " ";
				expText << setprecision(2);
				expText << exp << endl;
				totalExposure += exp;
				}
			else {
				status = 1005;
				break;
				}
			}
		else
			cout << "No intervals selected" << endl;
		beginTime = endTime;
		endTime += deltaT;
		if (tmax<endTime)
			endTime = tmax;
		} while (beginTime<tmax);
	expText.close();
	fits_close_file(localFits, &status);
	cout << "Total Exposure: " << totalExposure << endl;
	}
remove(s_selectionFileName);
remove(s_localFileName);

if (status) {
	if (status != 105)
		remove(outfile);
	cout << "AG_expmapgenT.................... exiting ERROR: " << status;
	if (status>0)
		fits_report_error(stdout, status);
	cout << endl;
	}
else
	cout << "AG_expmapgenT.................... exiting SUCCESS" << endl;

return status;
}
