////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       excalibur
//       Release: V0.0 -  8/Dec/2005
//       Contributors: 
//       Author: Andrea Bulgarelli (IASF_Bologna), derivated 
//	 from AG_expmapgen (Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano))
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       28/Sept/2007
//                      First release: V1.0
//       		Author: Andrea Bulgarelli (IASF-Bologna)
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


int AG_intmapgen(char * expfile, char *  outfile, char * ctsfile) {
	

	int status = 0;
	double x0 = 0, y0 = 0;	
	long outpixel[2];
	int bitpix   =  DOUBLE_IMG;
	int naxis;   
	long naxes[3] ;

	double l0 = 0, b0 = 0;
	double dl = 0, db = 0;	
	
	long nx = 0, ny = 0;

	AgileMap ExpMap(ctsfile);
	fitsfile * expFits;
	fitsfile * ctsFits;
	fitsfile * outFits;
	
	if ( fits_open_file(&expFits, expfile, READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n",  expfile);
		return status;
		}
	if ( fits_open_image(&ctsFits, ctsfile, READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", ctsfile);
		return status;
		}		
		
			
	if ( fits_create_file(&outFits, outfile, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", outfile);
		return status;
		}	

	fits_copy_file(expFits, outFits, 1, 1, 1, &status);
	
	fits_get_img_param(ctsFits, 3, &bitpix, &naxis, naxes, &status);
	fits_read_key(ctsFits,TLONG,"NAXIS1",&nx,NULL,&status);
	fits_read_key(ctsFits,TLONG,"NAXIS2",&ny,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CRVAL1",&l0,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CDELT1",&dl,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CRPIX1",&x0,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CRVAL2",&b0,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CDELT2",&db,NULL,&status);
	fits_read_key(ctsFits,TDOUBLE,"CRPIX2",&y0,NULL,&status);
	char * str7 =  "photons ** (cm**2 s sr)**(-1)";	
	fits_update_key(outFits, TSTRING,  "BUNIT", str7, NULL, &status);	
	
	long nrows = ExpMap.Rows();
	long ncols = ExpMap.Cols();
	
	long intensitypixel[3] = {1,1,1};

	/// const double conversion = dl * D2R * db * D2R;
	const double conversion = dl * DEG2RAD * db * DEG2RAD;


	for (outpixel[0]=1;outpixel[0]<=nrows;outpixel[0]++)
		for(outpixel[1]=1;outpixel[1]<=ncols;outpixel[1]++){
			double pixelsCts = 0.0;
			double pixelsExp = 0.0;
			double pixelsInt = 0.0;
			intensitypixel[0] = outpixel[0];
			intensitypixel[1] = outpixel[1];
			
			if (intensitypixel[1] > ny)
				intensitypixel[1] = ny;
			else if (intensitypixel[1] < 1)
				intensitypixel[1] = 1;
			
			pixelsExp = 0.0;
			fits_read_pix(expFits, TDOUBLE, intensitypixel, 1, NULL, &pixelsExp, NULL, &status) ;

			fits_read_pix(ctsFits, TDOUBLE, intensitypixel, 1, NULL, &pixelsCts, NULL, &status) ;
			
			if(pixelsExp != 0)
				pixelsInt =  pixelsCts / pixelsExp; 
			else
				pixelsInt = 0.0;
			
/*			if(pixelsInt > 0.01)				
				cout << intensitypixel[0] << ", " << intensitypixel[1] << " " << pixelsCts << " / " << pixelsExp << " = " << pixelsInt << endl;*/
			
			if(pixelsInt < 0 || pixelsInt > 0.012)
				pixelsInt = 0;
			fits_write_pix(outFits, TDOUBLE, outpixel, 1, &pixelsInt, &status) ; 
			if (status) throw;
		}

	fits_close_file(expFits, &status);
	fits_close_file(ctsFits, &status);
	fits_close_file(outFits, &status);

	return status;
}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char expfile[FLEN_FILENAME];
	char ctsfile[FLEN_FILENAME];	
	char outfile[FLEN_FILENAME];

	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("expfile", expfile);
	status = PILGetString("ctsfile", ctsfile);	
	status = PILGetString("outfile", outfile);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_intmapgen.cpp v.1.0 - 25/09/2007 - A.B. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Enter exposure file name = "<<expfile<< endl;
	cout << "Enter output file name = " <<outfile<< endl;
	cout << "Enter counts file name = "<< ctsfile << endl;
	cout << " "<< endl;
	cout << " "<< endl;	
	

	cout << "AG_intmapgen...............................starting"<< endl;		
	if (status == 0)	
		status = AG_intmapgen(expfile, outfile, ctsfile);
	cout << "AG_intmapgen............................... exiting"<< endl;		
	if (status) {
		if (status != 105) {
			if (outfile[0] == '!')
				remove(outfile+1);
			else
				remove(outfile);
			}
		printf("AG_intmapgen..................... exiting AG_intmapgen ERROR:");		
		fits_report_error(stdout, status);	
		return status;			
		}			
	else {
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_intmapgen........... exiting #################\n");
		printf("#################################################################\n\n\n");					
		}			
	
	return status;
}

