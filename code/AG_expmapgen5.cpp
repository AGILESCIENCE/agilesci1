////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_expmapgen5
//       Release: V0.1 - 2013-07-24
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
//         V0.1 -  2013-07-24
//         Author: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
//               Andrea Bulgarelli (IASF-Bologna), Tomaso Contessi (Nuove Idee sas)
//         Merging AG_expbatch e AG_exmapgen0
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
/// #include "AG_misclib.h"
/// #include "AG_calfiles.h"

#include "wcshdr.h"
#include "wcsmath.h"
#include "wcstrig.h"
#include "sph.h"

#include "pil.h"
#include "GenmapParams.h"
// #include "Rotation.h"
#include "Intervals.h"
#include "CalibUtils.h"
#include "MathUtils.h"


using namespace std;


inline double AG_expmapgen_area(double xbin, double ybin, double theta, ProjType projection)
{
    if (projection==ARC)
        return 0.0003046174197867085688996857673060958405 * xbin * ybin * Sinaa(0.0174532925199432954743716805978692718782 * theta);
    else
        return 0.0003046174197867085688996857673060958405 * xbin * ybin;
}


/** Unused at least for now
void SkyToInst(const Rotation& satRot, double l, double b, double& theta, double& phi)
{
double targetRa;
double targetDec;
eulerold(l, b, &targetRa, &targetDec, 2);
double cosRa = cos(targetRa*DegToRad);
double cosDec = cos(targetDec*DegToRad);
double sinRa = sin(targetRa*DegToRad);
double sinDec = sin(targetDec*DegToRad);
Vect target(cosDec*cosRa, sinRa*cosDec, sinDec);
target /= satRot;
theta = acos(target.y)*RadToDeg;
phi = atan2(target.z, -target.x)*RadToDeg;
}
*/

static void ReadFitsCol(fitsfile* file, double* array, const char* colName, long rowOffs, long rowCount, int* status)
{
    int colNum;
    fits_get_colnum(file, 1, const_cast<char*>(colName), &colNum, status);
    fits_read_col(file, TDOUBLE, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}

static void ReadFitsCol(fitsfile* file, short* array, const char* colName, long rowOffs, long rowCount, int* status)
{
    int colNum;
    fits_get_colnum(file, 1, const_cast<char*>(colName), &colNum, status);
    fits_read_col(file, TSHORT, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}


/// Some variable were accessed but not used in previous versions.
/// Their definition and usage has been put under ifdef __UNUSED__VARS__

/// #define __UNUSED__VARS__

#ifdef __UNUSED__VARS__
#define ALLOCATE_ARRAYS(count) \
	double* evttime = new double[count];\
	double* q1 = new double[count];\
	double* q2 = new double[count];\
	double* q3 = new double[count];\
	double* q4 = new double[count];\
	double* ra_y = new double[count];\
	double* dec_y = new double[count];\
	double* psi = new double[count];\
	double* gp = new double[count];\
	double* livetime = new double[count];\
	short*  mode = new short[count];\
	short*  log_status = new short[count];\
	double* earth_ra = new double[count];\
	double* earth_dec = new double[count];\
	short*  phase = new short[count]

#else

#define ALLOCATE_ARRAYS(count) \
	double* evttime = new double[count];\
	double* ra_y = new double[count];\
	double* dec_y = new double[count];\
	double* psi = new double[count];\
	double* gp = new double[count];\
	double* livetime = new double[count];\
	double* earth_ra = new double[count];\
	double* earth_dec = new double[count];\
	short*  phase = new short[count]

#endif


#ifdef __UNUSED__VARS__

#define DELETE_ARRAYS() \
	delete[] evttime;\
	delete[] q1;\
	delete[] q2;\
	delete[] q3;\
	delete[] q4;\
	delete[] ra_y;\
	delete[] dec_y;\
	delete[] psi;\
	delete[] gp;\
	delete[] livetime;\
	delete[] mode;\
	delete[] log_status;\
	delete[] earth_ra;\
	delete[] earth_dec;\
	delete[] phase

#else

#define DELETE_ARRAYS() \
	delete[] evttime;\
	delete[] ra_y;\
	delete[] dec_y;\
	delete[] psi;\
	delete[] gp;\
	delete[] livetime;\
	delete[] earth_ra;\
	delete[] earth_dec;\
	delete[] phase

#endif

int pil_chridx(char *p, char c)
{   int  l;

    for (l = 0; p[l]; l++)  if (c == p[l]) return(l);
    return(-1);
}


int	pil_curly_append(char **p, char *s, int nchars)
{   int	l, nl;
    char	*np;

    if (NULL == p)  return(PIL_NUL_PTR);
    if (nchars < 0) return(PIL_BAD_ARG);
    if (0 == nchars)
    {
        if (NULL != *p) return(PIL_OK);
    }
    else
    {
        if (NULL == s)  return(PIL_NUL_PTR);
    }

    if (NULL == *p)  {
        l = 0;
    }
    else  {
        l = strlen(*p);
    }

    nl = nchars + l + 1;						/* + 1 - for EOS */

    if (NULL == *p)  {
        np = (char *)PIL_malloc(nl);
    }
    else  {
        np = (char *)PIL_realloc((void *)*p, nl);
    }
    if (NULL == np)
    {   if (NULL != *p)  {
            PIL_free(*p);    /* on error deallocate string */
            *p = NULL;
        }
        return(PIL_NO_MEM);
    }

    if (nchars > 0) memcpy(np + l, s, nchars);			/* append new string */
    np[nl - 1] = 0;						/* signal EOS */

    *p = np;
    return(PIL_OK);
}


int pil_curly_expand(char *src, char **dst)
{   char		ev[PIL_CURLYSIZE];
    int		sl, l, l2, l3, tl, r;
    char		*p;


    if (NULL == src) return(PIL_NUL_PTR);
    if (NULL == dst) return(PIL_NUL_PTR);
    *dst = NULL;

    sl = tl = 0;
    for (;;)
    {   if (-1 == (l = pil_chridx(src + sl + tl, '$'))) break;	/* no '$' from current pos, so terminate */
        if ('{' != src[sl + tl + l + 1])
        {   tl += l + 1;
            continue;						/* '$' not followed by '{' */
        }

        l += tl;							/* total number of bytes _BEFORE_ '${...}' */
        if (-1 == (l2 = pil_chridx(src + sl + l + 2, '}'))) break;  /* no closing '}' */
        if (PIL_OK != (r = pil_curly_append(dst, src + sl, l))) return(r);  /* copy part _BEFORE_ '$' */
        sl += l + 2;						/* position on 1st char after '${' */

        if (l2 > 0)						/* if non-empty env.var name */
        {   memcpy(ev, src + sl, l2);				/* then copy its name to the temporary buffer */
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

inline void addexpval(
    long i, long ii,
    long mxdim,
    const vector<int>& aitstatus,
    const vector<double>& lng,
    const vector<double>& lat,
    double lp, double bp,
    ExpGenParams& params,
    double learth, double bearth,
    double time,
    double* A,
    AeffGridAverage* raeffArr,
    const vector<double>& area)
{
    long element = i * mxdim + ii;
    if (aitstatus[element] == 0 && params.AlbTest(lng[element], lat[element], learth, bearth)) {
        double theta = SphDistDeg(lng[element], lat[element], lp, bp);
        double phi = 0.0;
        long numout = params.maps.size();
        for (long ra = 0; ra < numout; ra++)
            if (params.FovTest(ra, theta))
                A[ra*mxdim*mxdim+element] += 1e-3*time*(raeffArr[ra].AvgVal(theta, phi))*area[element];
    }
}

/**
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
*/

static string LogExprString(const Interval& intv, int phasecode, int timeStep)
{
    stringstream str(ios_base::out);
    str.precision(6);
    str << "TIME >= " << fixed << intv.Start() << " && TIME < " << intv.Stop();
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



int excalibur(ExpGenParams& params)
{
    char buffer[1024];
    int status = 0;
    int binstep = params.binstep;
    long long numout = params.maps.size();

    double lp = 0, bp = 0;
    double learth, bearth;
    double lp0 = 0, bp0 = 0, gp0 = 0;
    double x = 0, y = 0;
    long i = 0, ii = 0;
    long mxdim= params.mxdim; // dimension (in pixels) of the map
    double theta, phi = 0;
    double eul[5];

    int intvCount = params.intervals.Count();
    double** matArr = new double*[intvCount];

    long numsquare = mxdim * mxdim;
    for (int i=0; i<intvCount; ++i) {
        matArr[i] = new double[numout*numsquare];
        double* A = matArr[i];
        for (long long ra = 0; ra < numout*numsquare; ra++)
            A[ra] = 0.0;
    }
    /** Replaced by above
    	double A[numout*numsquare];
    	for (long long ra = 0; ra < numout*numsquare; ra++)
    		A[ra] = 0.0;
    */

    vector<double> lng(numsquare);
    double lng0, lat0;
    vector<double> lat(numsquare);
    vector<double> area(numsquare);
    vector<int> aitstatus(numsquare);
    int ait0;
    struct prjprm *prj = new(prjprm);
    prjini(prj);
    switch (params.projection) {
    case ARC:
        for (i = 0; i <= mxdim-1; ++i) {
            x = -(params.mdim/2)+(params.mres*(i+0.5));
            for (ii = 0; ii <= mxdim-1; ++ii) {
                aitstatus[i*mxdim+ii] = 0;
                y = -(params.mdim/2)+(params.mres*(ii+0.5));
                theta = 90-sqrt(x*x+y*y);
                phi = atan2d(-y, -x);
                eul[0] = params.la;
                eul[1] = 90-params.ba;
                eul[2] = params.lonpole;
                eul[3] = cosd(eul[1]);
                eul[4] = sind(eul[1]);
                sphx2s(eul, 1, 1, 0, 0, &phi, &theta, &lng0, &lat0);
                lng[i*mxdim+ii] = lng0;
                lat[i*mxdim+ii] = lat0;
                area[i*mxdim+ii] = AG_expmapgen_area(params.mres, params.mres, 90-theta, params.projection);
            }
        }
        break;
    case AIT:
        prj->r0 = 180/M_PI;
        prj->phi0 = params.la;
        prj->theta0 = params.ba;
        aitset(prj);
        for (i = 0; i <= mxdim-1; ++i) {
            y = -(params.mdim/2)+(params.mres*(i+0.5));
            for (ii = 0; ii <= mxdim-1; ++ii) {
                x = (params.mdim/2)-(params.mres*(ii+0.5));
                aitx2s(prj, 1, 1, 0, 0, &x, &y, &lng0, &lat0, &ait0);
                lng[i*mxdim+ii] = lng0;
                lat[i*mxdim+ii] = lat0;
                aitstatus[i*mxdim+ii] = ait0;
                area[i*mxdim+ii] = AG_expmapgen_area(params.mres, params.mres, 0, params.projection);
            }
        }
        break;
    }
    int bitpix = DOUBLE_IMG; /* 16-bit unsigned short pixel values */
    long naxis = 2;  /* 2-dimensional image */
    long naxes[2] = { mxdim, mxdim };
    long n = 0;
    double time = 0;

    FILE *fp = fopen(params.logfile, "r");
    if (!fp) {
        cerr << "Error reading " << params.logfile << endl;
        return 104;
    }

    /**
    	char expr[1024];
    	strcpy(expr,params.logexpr().c_str());
    	char expr2[4096];
    	char timestepstring[128];
    	sprintf(timestepstring, "%d", int(params.timestep)); /// zzz this is a double

    	sprintf(expr2,"%s && ((#ROW == 1) || (#ROW == (#ROW/%s) *%s))", expr, timestepstring, timestepstring);
    	cout << expr2 << endl;
    */


    int hdutype = 0;
    int  find = 0;
    char *name=new char[FLEN_FILENAME];
    char *nname=new char[FLEN_FILENAME];
    char *tempfilenamegz=new char[FLEN_FILENAME];
    char *tempfilename=new char[FLEN_FILENAME];
    char *command=new char[FLEN_FILENAME];
    char *fileext=new char[FLEN_FILENAME];

    double t1 = 0, t2 = 0;

    /// AeffGridAverage *raeff[20]; /// zzz allocate dynamically
    AeffGridAverage *raeffArr = new AeffGridAverage[numout];

    bool hasEdp = strcmp(params.edpFileName, "None");
    for (long ra = 0 ; ra < numout ; ra++) {
        int status = raeffArr[ra].Read(params.sarFileName, params.maps[ra].emin, params.maps[ra].emax, params.maps[ra].index);
        if (status) {
            cerr << "Error reading " << params.sarFileName << endl;
            return status;
        }
        /// AeffGridAverage* raeffPtr = new AeffGridAverage(params.sarFileName, params.maps[ra].emin, params.maps[ra].emax, params.maps[ra].index);
        if (hasEdp)
            /// raeffPtr->LoadEdp(params.edpFileName);
            status = raeffArr[ra].LoadEdp(params.edpFileName);
        /// raeff[ra] = raeffPtr;
        if (status) {
            cerr << "Error reading " << params.edpFileName << endl;
            return status;
        }
    }

    while (fgets(buffer , 40960, fp)) {
        sscanf(buffer, "%s %lf %lf", name, &t1, &t2);
        Interval fileIntv(t1, t2);
        Intervals fileIntersect = Intersection(params.intervals, fileIntv);
        if (fileIntersect.Count()) {

            if (fileIntersect.Count()>1)
                cout << "Intervals: " << fileIntersect.Count() << endl;


            /// if ( ((t2 > params.tmin && t2 < params.tmax)  || (t1 > params.tmin && t1 < params.tmax)  || (t1 <= params.tmin && t2 >= params.tmax))) {
            pil_curly_expand(name, &nname);
//			cout<<nname<<endl;
            tmpnam(tempfilename);
            for (int ffi=0; ffi<20; ffi++)
                fileext[ffi]='A'+rand()%26;
            fileext[20]=0;
            strcat(tempfilename,fileext);
            strcpy(tempfilenamegz, tempfilename);
            strcat(tempfilenamegz,".gz");
            cout  << nname << endl;
            std::ifstream ifs(nname, std::ios::binary);
            std::ofstream ofs(tempfilenamegz, std::ios::binary);
            ofs << ifs.rdbuf();
            ifs.close();
            ofs.close();

            strcpy(command, "gunzip ");
            strcat(command, tempfilenamegz);
            cout << "Command: " << command << endl;
            system(command);
//			cout << "Closing " << nname << endl;
//			fits_close_file(tempFits2, &status);
//			cout << "Opening file " << tempfilename << endl;

            for (int intvIndex = 0; intvIndex<fileIntersect.Count(); ++intvIndex) {
                Interval thisIntv(fileIntersect[intvIndex]);
                int intvGlobIndex = params.intervals.IndexOf(thisIntv.Start());
                if (intvGlobIndex<0)
                    continue; /// This should never happen
                double* A = matArr[intvGlobIndex]; /// Select the matrix for this time interval

                fitsfile *tempFits;
                if ( fits_open_file(&tempFits, tempfilename, READWRITE, &status) != 0 ) {
                    cerr << "Error opening file " << name << endl;
                    return status;
                }
//			cout << "Opened file " << tempfilename << endl;
                fits_movabs_hdu(tempFits, 2, &hdutype, &status);
//			cout << "Selecting rows from " << tempfilename << endl;
                string selExpr = LogExprString(thisIntv, params.phasecode, params.timestep);
                {
                    char tmpSel[1024];
                    strcpy(tmpSel, selExpr.c_str());
                    fits_select_rows(tempFits, tempFits, tmpSel, &status);
                    //			cout << "Rows from " << tempfilename << " selected" << endl;
                }
                find++;

//			fits_movabs_hdu(tempFits, 2, NULL, &status);
                long allnrows;
                fits_get_num_rows(tempFits, &allnrows, &status);
//			cout << allnrows << " rows total" << endl;

                long rowblockincrement = 65536;
                fits_get_rowsize(tempFits, &rowblockincrement, &status);
//			cout << "Row block size = " << rowblockincrement << endl;
//			cout << "Allocating arrays" << endl;

                ALLOCATE_ARRAYS(rowblockincrement);

                long* change = new long[rowblockincrement];

                for (long rowblockzero=0 ; rowblockzero < allnrows ; rowblockzero += rowblockincrement) {

                    long nrows = (rowblockzero + rowblockincrement < allnrows) ? rowblockincrement : (allnrows - rowblockzero);

//	cout << "Reading " << nrows << "rows" << endl;
                    ReadFitsCol(tempFits, evttime, "TIME", rowblockzero, nrows, &status);
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
                        if ((phase[k-1] != phase[k]) || params.YTolTest(ra_y[k-1], dec_y[k-1], ra_y0, dec_y0) || params.RollTolTest(&psi[k-1]) || params.EarthTolTest(earth_ra[k-1], earth_dec[k-1], earth_ra0, earth_dec0)) {
                            change[count++] = k;
                            earth_ra0 = earth_ra[k-1];
                            earth_dec0 = earth_dec[k-1];
                            ra_y0 = ra_y[k-1];
                            dec_y0 = dec_y[k-1];
                        }
                    }
                    if (count == 0)
                        change[0] = nrows;
                    //	cout << endl << rowblockzero << " Count = " << count << endl;


                    long lowrow = 0, highrow = 0;

                    for (long k = 0; k<=count; ++k) {
                        time = 0;
                        if ( k == 0 ) {
                            lowrow = 0;
                            highrow = change[0];
                        }
                        else if ( k == count ) {
                            lowrow = change[k-1];
                            highrow = nrows;
                        }
                        else {
                            lowrow = change[k-1];
                            highrow = change[k];
                        }

                        for (n = lowrow; n<highrow; n++)
                            time += (livetime[n] * params.timestep);

                        //		euler(ra_y[lowrow], dec_y[lowrow], psi[lowrow], &lp, &bp, gp+lowrow);

                        /// eulerold(ra_y[lowrow], dec_y[lowrow], &lp, &bp, 1);
                        Euler(ra_y[lowrow], dec_y[lowrow], &lp, &bp, 1);

                        /// eulerold(earth_ra[lowrow], earth_dec[lowrow], &learth, &bearth, 1);
                        Euler(earth_ra[lowrow], earth_dec[lowrow], &learth, &bearth, 1);

                        if (k==0 && rowblockzero==0) {
                            lp0 = lp;
                            bp0 = bp;
                            // gp0=gp[lowrow]*R2D;
                        }

                        // reportfile << "LP = " << lp << " BP = " << bp << " LonGPole= " << R2D*gp[lowrow]<< endl;
                        // reportfile << "Earth (theta, phi) = " << learth << ", " << bearth << endl;
                        // cout << binstep << " " << mxdim << endl;

                        /// Rotation satRot(q4[lowrow], q1[lowrow], q2[lowrow], q3[lowrow]);

                        // double factor;
                        for (i = 0; i <= mxdim-2; i+= binstep) {
                            for (ii = 0; ii <= mxdim-2; ii+=binstep)
                                addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, params, learth, bearth, time, A, raeffArr, area);
                            ii = mxdim - 1;
                            addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, params, learth, bearth, time, A, raeffArr, area);
                        }
                        i = mxdim - 1;
                        for (ii = 0; ii <= mxdim-2; ii+=binstep)
                            addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, params, learth, bearth, time, A, raeffArr, area);
                        ii = mxdim - 1;
                        addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, params, learth, bearth, time, A, raeffArr, area);
                    }
                }
                //	cout << "Deleting " << tempfilename << endl;
                fits_delete_file(tempFits, &status);
//			cout << "Deleting arrays" << endl;
                DELETE_ARRAYS();
//			cout << "Arrays deleted" << endl;
                if (status)
                    return status;

            }
        }
    }

    /**
    for (long i = 0; i < numout ; i++)
    	delete raeff[i];
    */
    delete[] raeffArr;

    if (find == 0)
        return 1005;

    fclose(fp);
    cout << "Log file closed" << endl;
    delete [] name;
    delete [] nname;
    delete [] tempfilenamegz;
    delete [] tempfilename;
    delete [] command;
    delete [] fileext;


    for (int matIndex=0; matIndex<intvCount; ++matIndex) {
        double* A = matArr[matIndex]; /// Select the matrix for a time interval
        Interval intv = params.intervals[matIndex]; /// the time interval

// 	#ifdef AAAA
        long nrows = mxdim;
        long ncols = mxdim;
        long outpixel[2];
        double pixel = 0.0;
        //AB
        //	1	3
        //	4	2
        double pixel2 = 0.0;
        double pixel3 = 0.0;
        double pixel4 = 0.0;

        cout << "Beginning interpolation" << endl;
        if(binstep > 1) {
            for (long long ra = 0; ra < numout ; ra++) {
                int step0 = binstep;
                bool end0 = true;
//				for (outpixel[0]=nrows-binstep;end0;outpixel[0]-=binstep) {
                for (outpixel[0]=1; end0; outpixel[0]+=binstep) {
                    bool end1 = true;
                    int step1 = binstep;
// 		  if (outpixel[0] < 1) {
// 		    step0 = binstep + outpixel[0]-1;
// 		    outpixel[0] = 1;
// 		    end0 = false;
                    if (outpixel[0] + binstep > nrows) {
                        step0 = nrows - outpixel[0];
                        end0 = false;
                    }
                    for(outpixel[1]=1; end1; outpixel[1]+=binstep) {
                        if (outpixel[1] + binstep > ncols) {
                            step1 = ncols - outpixel[1];
                            end1 = false;
                        }
// 			  cout << outpixel[1]-1 << "," << outpixel[0]-1 << ":" << end1 << "," << end0 << ":" << step1 << "," << step0 << ";";
                        pixel = A[ra*numsquare+(outpixel[1]-1)*mxdim+(outpixel[0]-1)];
                        if(pixel != 0.0 && step0 > 0 && step1 > 0) {
                            pixel2 = A[ra*numsquare+(outpixel[1]-1 + step1)*mxdim+(outpixel[0]-1 + step0)];
                            pixel3 = A[ra*numsquare+(outpixel[1]-1 + step1)*mxdim+(outpixel[0]-1)];
                            //si lavora sul primo triangolo
                            long y1 = outpixel[0];
                            long x1 = outpixel[1];
                            double z1 = pixel;
                            long y2 = outpixel[0] + step0;
                            long x2 = outpixel[1] + step1;
                            double z2 = pixel2;
                            long y3 = outpixel[0];
                            long x3 = outpixel[1] + step1;
                            double z3 = pixel3;
                            double a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
                            double b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
                            long c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
                            // 				cout << "1 " << x1 << " " << y1 << " " << z1 << endl;
                            // 				cout << "2 " << x2 << " " << y2 << " " << z2 << endl;
                            // 				cout << "3 " << x3 << " " << y3 << " " << z3 << endl;
                            for(int i=y1; i<=y2; i++)
                                for(int j=x1; j<=x2; j++) {
                                    if(i-y1<=j-x1) {
                                        double p = (c * z1 - a * (j - x1) - b * (i - y1)) / c;
                                        long outpixel2[2];
                                        long x, y;
                                        outpixel2[0] = y = i;
                                        outpixel2[1] = x = j;
                                        if (x < 1 || x > mxdim || y < 1 || y > mxdim )
                                            cout << "x " << x << " " << y << " " << p << endl;
                                        A[ra*numsquare+(outpixel2[1]-1)*mxdim+(outpixel2[0]-1)] = p;
                                    }
                                }

                            pixel4 = A[ra*numsquare+(outpixel[1] - 1)*mxdim+(outpixel[0]-1 + step0)];
//					cout << pixel << "," <<  pixel2 << "," <<  pixel3 << "," <<  pixel4 << " ; ";
                            //si lavora sul secondo triangolo
                            y3 = outpixel[0] + step0;
                            x3 = outpixel[1];
                            z3 = pixel4;
                            a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
                            b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
                            c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
                            // 				cout << "4 " << x3 << " " << y3 << " " << z3 << endl;
                            for(int i=y1; i<=y2; i++)
                                for(int j=x1; j<=x2; j++) {
                                    if(i-y1>j-x1) {
                                        double p = (c * z1 - a * (j - x1) - b * (i - y1)) / c;
                                        long outpixel2[2];
                                        long x, y;
                                        outpixel2[0] = y = i;
                                        outpixel2[1] = x = j;
                                        if (x < 1 || x > mxdim || y < 1 || y > mxdim )
                                            cout << "x " << x << " " << y << " " << p << endl;
                                        A[ra*numsquare+(outpixel2[1]-1)*mxdim+(outpixel2[0]-1)] = p;
                                    }
                                }
                        }
                    }
                }
            }
        } //endif
        cout << "Ending interpolation" << endl;
        fitsfile * mapFits;
        for (long long i=0; i<numout; i++) {
            double obtMin = int(intv.Start());
            double obtMax = int(intv.Stop());
            char outfile[2048] = "!";
			if(numout > 1) {
				if (params.outfile[0]=='!')
					sprintf(outfile+1, "TA%d_TB%d_E%d-%d_TH%d-%d_SI-%g_%s",
							int(obtMin), int(obtMax),
							int(params.maps[i].emin), int(params.maps[i].emax),
							int(params.maps[i].fovradmin), int(params.maps[i].fovradmax),
							params.maps[i].index, params.outfile+1);
				else
					sprintf(outfile, "TA%d_TB%d_E%d-%d_TH%d-%d_SI-%g_%s",
							int(obtMin), int(obtMax),
							int(params.maps[i].emin), int(params.maps[i].emax),
							int(params.maps[i].fovradmin), int(params.maps[i].fovradmax),
							params.maps[i].index, params.outfile);
			} else {
				if (params.outfile[0]=='!')
					sprintf(outfile+1, "%s", params.outfile+1);
				else
					sprintf(outfile, "%s", params.outfile);
			}
            /// string outfile = params.maps[i].outfile;

            cout << "Creating file " << outfile << endl;
            if ( fits_create_file(&mapFits, outfile, &status) != 0 ) {
                cerr << "Error opening file " << outfile << endl;
                return status;
            }
            cout << "Created file " << outfile << endl;

            cout << "Creating image in " << outfile << endl;
            fits_create_img(mapFits, bitpix, naxis, naxes, &status);
//  int fits_write_img_[byt, sht, usht, int, uint, lng, ulng, lnglng, flt, dbl] /      ffppr[b,i,ui,k,uk,j,uj,jj,e,d]      (fitsfile *fptr, long group, LONGLONG firstelem, LONGLONG nelements,       DTYPE *array, > int *status);
//int fits_write_pix / ffppx(fitsfile *fptr, int datatype, long *fpixel, LONGLONG nelements,DTYPE *array, int *status);
            long startpixel[2] = {1,1};
            fits_write_pix(mapFits, TDOUBLE, startpixel, numsquare, &A[i * numsquare], &status);
//			 fits_write_img_dbl(mapFits, 0, 1, mxdim * mxdim, &A[i * mxdim * mxdim], &status);
            params.write_fits_header(i, mapFits, status);

            fits_update_key(mapFits, TDOUBLE, "TSTART", &obtMin, "[OBT]first event time", &status);
            fits_update_key(mapFits, TDOUBLE, "TSTOP", &obtMax, "[OBT]last event time", &status);

            char timeSys[] = "TT";
            char timeUnit[] = "s";
            double zero = 0;
            fits_update_key(mapFits, TSTRING, "TIMESYS", timeSys, "TT", &status);
            fits_update_key(mapFits, TSTRING, "TIMEUNIT", timeUnit, "TT", &status);
            fits_update_key(mapFits, TDOUBLE, "TZERO", &zero, timeSys, &status);

            fits_update_key(mapFits, TDOUBLE, "SC-Z-LII", &lp0, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "SC-Z-BII", &bp0, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "SC-LONPL", &gp0, NULL, &status);
            cout << "Closing " << outfile << endl;
            fits_close_file(mapFits, &status);
        }
    }

    for (int matIndex=0; matIndex<intvCount; ++matIndex)
        delete[] matArr[matIndex];
    delete[] matArr;

//	for (long ra=0; i<numout; i++) {
//		for (long i=0; i<mxdim; i++)
//			delete A[ra][i];
//		delete A[ra];
//	}
//	delete A;
    return status;
}



class AppScreen
{
public:
    AppScreen()
    {
        cout << "#################################################################" << endl;
        cout << "#### AG_expmapgen5 - A.C., T.C., A.T. A.B.                  #####" << endl;
        cout << "#################################################################" << endl;
        cout << endl;
    }

    ~AppScreen()
    {
        cout << endl;
        cout << "#################################################################" << endl;
        cout << "##########  Task AG_expmapgen5........... exiting ###############" << endl;
        cout << "#################################################################" << endl;
    }
};



int main(int argc, char* argv[])
{
    AppScreen appScreen;

    ExpGenParams params;
    if (!params.Load(argc, argv))
        return -1;

    int status = excalibur(params);

    /*
    	int status = 0;
    	int numpar=0;

    	AG_batchparams params(AG_batchparams::EXP, argc, argv, numpar, status);

    	cout << "AG_expmapgen...............................starting"<< endl;
    	if (status == 0)
    		status = excalibur(params);
    	cout << "AG_expmapgen............................... exiting"<< endl;
    */
    if (status) {
        /**
        if (status != 105)
        	for (unsigned long i = 0; i < params.maps.size() ; i++ )
        		remove(params.maps[i].outfile.c_str());
        */
        printf("AG_expmapgen..................... exiting AG_expmapgen ERROR:");
        fits_report_error(stdout, status);
        return status;
    }
    return status;
}
