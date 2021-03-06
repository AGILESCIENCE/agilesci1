////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Norm
//       Authors: Andrea Zoli (INAF/IASF Bologna)
//
// INPUT
//       A map
//
// OUTPUT
//       The input map normalized
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


#include <fitsio.h>
#include <pil.h>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

const char* startString = {
"###################################################\n"
"###      AG_norm B25 v1.2.0 - A.Z               ###\n"
"###################################################\n"
};

const char* endString = {
"###################################################\n"
"###      AG_norm B25 ended successfully ###########\n"
"###################################################\n"
};
const char* errString = {
"###################################################\n"
"###     AG_norm B25 ended with errors #############\n"
"###################################################\n"
};


void exitOnFitsIOError(int status) {
    if (status) {
        cerr << endl << "FITS ERROR!";
        fits_report_error(stderr, status);
        cerr << endl << errString << endl;
        exit(status);
    }
}

void exitOnPilError(int status) {
    if (status != PIL_OK) {
        cerr << endl << "PIL ERROR!" << endl;
        cerr << "Error code: " << status << endl;
        cerr << endl << errString << endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[]) {
    cout << startString << endl;
	cout << "Using " << *pilversion << endl;

    int status, numpar;
	char mapFilename[FLEN_FILENAME];
	char mapnormFilename[FLEN_FILENAME];

	// parse .par or command line options with PIL
	status = PILInit(argc,argv);
	exitOnPilError(status);
    status = PILGetNumParameters(&numpar);
	exitOnPilError(status);
	status = PILGetString("map", mapFilename);
	exitOnPilError(status);
	status = PILGetString("mapnorm", mapnormFilename);
	exitOnPilError(status);
	status = PILClose(status);
	exitOnPilError(status);

	cout << endl << "INPUT PARAMETERS:" << endl;
	cout << "Input map filename = " << mapFilename << endl;
	cout << "Normalized map filename = " << mapnormFilename << endl;

    fitsfile *ifd;
    fits_open_file(&ifd, mapFilename, READONLY, &status);
    exitOnFitsIOError(status);

    fitsfile *ofd;
    fits_create_file(&ofd, mapnormFilename, &status);
    exitOnFitsIOError(status);
    fits_copy_header(ifd, ofd, &status);
    exitOnFitsIOError(status);

    int numaxis;
    fits_get_img_dim(ifd, &numaxis, &status);
    exitOnFitsIOError(status);

    if(numaxis != 2) {
        cerr << "Error: only images with 2 dimensions are supported!" << endl;
        return EXIT_FAILURE;
    }

    long dimaxis[2];
    fits_get_img_size(ifd, 2, dimaxis, &status);
    exitOnFitsIOError(status);

    long npixels = dimaxis[0] * dimaxis[1];
    long pix1[3] = {1, 1, 1};

    double *inbuff = new double[npixels];
    fits_read_pix(ifd, TDOUBLE, pix1, npixels, NULL, inbuff, NULL, &status);
    exitOnFitsIOError(status);
    fits_close_file(ifd, &status);
    exitOnFitsIOError(status);

    // compute the sum
    double sum = 0.0;
    for (int y=0; y<dimaxis[0]; ++y) {
        for (int x=0; x<dimaxis[1]; ++x) {
            int idx = y * dimaxis[1] + x;
            sum += inbuff[idx];
        }
    }
    // normalize
    double *normbuff = new double[npixels];
    for (int y=0; y<dimaxis[0]; ++y) {
        for (int x=0; x<dimaxis[1]; ++x) {
            int idx = y * dimaxis[1] + x;
            normbuff[idx] = inbuff[idx] / sum;
        }
    }

    delete[] inbuff;
    fits_write_pix(ofd, TDOUBLE, pix1, npixels, normbuff, &status);
    delete[] normbuff;
    exitOnFitsIOError(status);

    fits_close_file(ofd, &status);
    exitOnFitsIOError(status);

    cout << endString << endl;

    return EXIT_SUCCESS;
}
