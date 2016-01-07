////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       Release: V0.0 - 07/Jan/2016
//       Contributors:
//       Author: Andrea Bulgarelli (IASF_Bologna), derivated from
//       AG_difmapgen (Andrea Bulgarelli (IASF_Bologna), Andrew Chen,
//       Alberto Pellizzoni, Alessio Trois (IASF-Milano))
//
// INPUT
//       Two 2d image as fits with the same number of axes
//
// OUTPUT
//       The difference between the two images
//
// FILE HISTORY
//       07/Jan/2016
//       First release: V1.0
//       Author: Andrea Zoli (IASF-Bologna)
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>

#include "fitsio.h"
#include "pil.h"

#include <AgileMap.h>
#include <MathUtils.h>

using namespace std;

int AG_difmapgen(char *map1file, char *map2file, char *outfile) {
	int status = 0;
	double x0 = 0, y0 = 0;
	long outpixel[2];
	int bitpix = DOUBLE_IMG;
	int naxis;
	long naxes[3];

	double l0 = 0, b0 = 0;
	double dl = 0, db = 0;

	long nx = 0, ny = 0;

	fitsfile *map1Fits;
	fitsfile *map2Fits;
	fitsfile *outFits;

	if (fits_open_file(&map1Fits, map1file, READONLY, &status) != 0) {
		printf("Errore in apertura file '%s'\n",  map1file);
		return status;
	}
	if (fits_open_image(&map2Fits, map2file, READONLY, &status) != 0) {
		printf("Errore in apertura file '%s'\n", map2file);
		return status;
	}
	if (fits_create_file(&outFits, outfile, &status) != 0) {
		printf("Errore in apertura file '%s'\n", outfile);
		return status;
	}

	fits_copy_file(map1Fits, outFits, 1, 1, 1, &status);

	fits_get_img_param(map2Fits, 3, &bitpix, &naxis, naxes, &status);
	fits_read_key(map2Fits,TLONG,"NAXIS1",&nx,NULL,&status);
	fits_read_key(map2Fits,TLONG,"NAXIS2",&ny,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CRVAL1",&l0,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CDELT1",&dl,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CRPIX1",&x0,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CRVAL2",&b0,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CDELT2",&db,NULL,&status);
	fits_read_key(map2Fits,TDOUBLE,"CRPIX2",&y0,NULL,&status);

	AgileMap Map(map2file);
	long nrows = Map.Rows();
	long ncols = Map.Cols();

	long pixel[3] = {1,1,1};

	for (outpixel[0]=1;outpixel[0]<=nrows;outpixel[0]++) {
		for (outpixel[1]=1;outpixel[1]<=ncols;outpixel[1]++) {
			double pixelsMap2 = 0.0;
			double pixelsMap1 = 0.0;
			double pixelsInt = 0.0;
			pixel[0] = outpixel[0];
			pixel[1] = outpixel[1];

			if (pixel[1] > ny)
				pixel[1] = ny;
			else if (pixel[1] < 1)
				pixel[1] = 1;

			fits_read_pix(map1Fits, TDOUBLE, pixel, 1, NULL, &pixelsMap1, NULL, &status);
			fits_read_pix(map2Fits, TDOUBLE, pixel, 1, NULL, &pixelsMap2, NULL, &status);

			pixelsInt =  pixelsMap2 - pixelsMap1;
/*
			if(pixelsInt > 0.01)
				cout << pixel[0] << ", " << pixel[1] << " " << pixelsMap2 << " - " << pixelsMap1 << " = " << pixelsInt << endl;
*/
			fits_write_pix(outFits, TDOUBLE, outpixel, 1, &pixelsInt, &status);
			if (status) throw;
		}
	}

	fits_close_file(map1Fits, &status);
	fits_close_file(map2Fits, &status);
	fits_close_file(outFits, &status);

	return status;
}


int main(int argc,char **argv)
{
	int status = 0, numpar = 0;
	char map1file[FLEN_FILENAME];
	char map2file[FLEN_FILENAME];
	char outfile[FLEN_FILENAME];

	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("map1file", map1file);
	status = PILGetString("map2file", map2file);
	status = PILGetString("outfile", outfile);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;
	cout << "#################################################################"<< endl;
	cout << "########## AG_difmapgen.cpp v.1.0 - 07/01/2016          #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Enter first map filename = "  << map1file << endl;
	cout << "Enter second map filename = " << map2file << endl;
	cout << "Enter output map filename = " << outfile << endl;
	cout << " "<< endl;
	cout << " "<< endl;

	cout << "AG_difmapgen...............................starting"<< endl;
	if (status == 0)
		status = AG_difmapgen(map1file, map2file, outfile);
	cout << "AG_difmapgen............................... exiting"<< endl;
	if (status) {
		if (status != 105) {
			if (outfile[0] == '!')
				remove(outfile+1);
			else
				remove(outfile);
		}
		printf("AG_difmapgen..................... exiting AG_difmapgen ERROR:");
		fits_report_error(stdout, status);
		return status;
	}
	else {
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_difmapgen........... exiting #################\n");
		printf("#################################################################\n\n\n");
	}

	return status;
}
