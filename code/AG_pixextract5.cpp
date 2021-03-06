////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG pixextract
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       Copyright (C) 2005-2019 AGILE Team. All rights reserved.
/*
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <sstream>
#include "pil.h"
#include "fitsio.h"

using namespace std;


// From gridutilities
const double obtlimit = 104407200.0;

static int addfile(fitsfile* iFile, const char* fileList, char* expr, double tmin, double tmax)
{
char buffer[40960];
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


int pixextract(
	const char* infile,
	char* logfile,
	const char* outfile,
	double tmin,
	double tmax,
	double tdur,
	double radius)
{
double l, b;
int i, status=0;
long nrows;
long minnrows = (long)tdur/16;
ifstream input(infile);
ofstream output(outfile);
while (!input.eof()) {
	input >> i >> b >> l;
	fitsfile* tempFits;
    stringstream ss;
    ss << "/tmp/AG_pixelextract5_" << getpid();
    std::string tmpstr = ss.str();
    char tempname[FLEN_FILENAME];
    strncpy(tempname, tmpstr.c_str(), FLEN_FILENAME-1);
	if (fits_create_file(&tempFits, tempname, &status) != 0) {
		printf("Errore in apertura file %s\n",tempname);
		return status;
		}
	char expr[1024];
	sprintf(expr, "TIME >= %f && TIME <= %f && LIVETIME > 0 && LOG_STATUS == 0 && MODE == 2 && angsep(%f,%f,ATTITUDE_GLON_Y,ATTITUDE_GLAT_Y) < %f", tmin, tmax, l, b, radius);
	status = addfile(tempFits, logfile, expr, tmin, tmax);
	if (status)
		return status;
	fits_get_num_rows(tempFits, &nrows, &status);
	if (nrows >= minnrows)
		output << i << " " << b << " " << l << endl;
	fits_delete_file(tempFits, &status);
	}
output.close();
input.close();
return status;
}


int main(int argc, char* argv[])
{
int status = 0;
int numpar=0;
char infile[FLEN_FILENAME];
char logfile[FLEN_FILENAME];
char outfile[FLEN_FILENAME];
double tmin, tmax, radius, tdur;

status = PILInit(argc, argv);
status = PILGetNumParameters(&numpar);
status = PILGetString("infile", infile);
status = PILGetString("logfile", logfile);
status = PILGetString("outfile", outfile);
status = PILGetReal("tmin", &tmin);
status = PILGetReal("tmax", &tmax);
status = PILGetReal("tdur", &tdur);
status = PILGetReal("radius", &radius);

if (status==0)
	status = pixextract(infile, logfile, outfile, tmin, tmax, tdur, radius);
if (status!=0 && status!=105)
	fits_report_error(stdout, status);
return status;
}

