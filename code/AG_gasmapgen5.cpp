////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_gasmapgen2
//       Release: V2.0 -  10/Oct/2010
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
//       10/Oct/2010
//                      Second release: V2.0
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstring>

#include <iostream>

#include "fitsio.h"
#include "pil.h"

#include <AgileMap.h>

using namespace std;


int AG_gasmapgen(char * expfile, char *  outfile, char * loresdiffusefile, char * hiresdiffusefile) {
	

//  	ofstream outtextstream("gasmapgen.out");
//	char errormessage[FLEN_ERRMSG];
//	char errorstatus[FLEN_STATUS];
	int status = 0;

	long outpixel[2];
	int bitpix   =  DOUBLE_IMG;
	double pixels = 0.0;
	double nulval = 1.0/0.0;
	int anynul = 0;

	double lox0 = 0, loy0 = 0;	
	int lonaxis;   
	long lonaxes[3] ;
	double lol0 = 0, lob0 = 0;
	double lodl = 0, lodb = 0;	
	long lonx = 0, lony = 0;


	int hiresstatus = 0;
	double hix0 = 0, hiy0 = 0;
	int hinaxis;   
	long hinaxes[3] ;
	double hil0 = 0, hib0 = 0;
	double hidl = 0, hidb = 0;	
	long hinx = 0, hiny = 0;

	AgileMap ExpMap(expfile);
	fitsfile * expFits;
	fitsfile * lodiffuseFits;
	fitsfile * hidiffuseFits;
	
	fitsfile * outFits;
	
	if ( fits_open_file(&expFits, expfile, READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n",  expfile);
		return status;
		}
	if ( fits_open_image(&lodiffuseFits, loresdiffusefile, READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", loresdiffusefile);
		return status;
		}		
	
	if ( fits_open_image(&hidiffuseFits, hiresdiffusefile, READONLY, &hiresstatus) != 0 ) {
		printf("Cannot read high res diffuse file '%s'\n", hiresdiffusefile);
	} else {
		fits_get_img_param(hidiffuseFits, 3, &bitpix, &hinaxis, hinaxes, &hiresstatus);
		fits_read_key(hidiffuseFits,TLONG,"NAXIS1",&hinx,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TLONG,"NAXIS2",&hiny,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CRVAL1",&hil0,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CDELT1",&hidl,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CRPIX1",&hix0,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CRVAL2",&hib0,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CDELT2",&hidb,NULL,&hiresstatus);
		fits_read_key(hidiffuseFits,TDOUBLE,"CRPIX2",&hiy0,NULL,&hiresstatus);
		if ( hiresstatus != 0 )
			printf("ERROR: '%s' will not be used because a FITS keyword is missing\n", hiresdiffusefile);
	}
		
	if ( fits_create_file(&outFits, outfile, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", outfile);
		return status;
		}	

	fits_copy_file(expFits, outFits, 1, 1, 1, &status);
	
	fits_get_img_param(lodiffuseFits, 3, &bitpix, &lonaxis, lonaxes, &status);
	fits_read_key(lodiffuseFits,TLONG,"NAXIS1",&lonx,NULL,&status);
	fits_read_key(lodiffuseFits,TLONG,"NAXIS2",&lony,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CRVAL1",&lol0,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CDELT1",&lodl,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CRPIX1",&lox0,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CRVAL2",&lob0,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CDELT2",&lodb,NULL,&status);
	fits_read_key(lodiffuseFits,TDOUBLE,"CRPIX2",&loy0,NULL,&status);

	char str7[] =  "(cm**2 s sr)**(-1)";	
	fits_update_key(outFits, TSTRING,  "BUNIT", str7, NULL, &status);	
	char keyword[FLEN_FILENAME];
	char comment[FLEN_FILENAME];
	strcpy(keyword, "GAS");
	fits_update_key(outFits, TSTRING,  "EXTNAME", keyword, NULL, &status);
	if (fits_read_key(lodiffuseFits,TSTRING,"DH_CONF_",keyword,comment,&status) == 0)
		fits_update_key(outFits,TSTRING,"DH_CONF_",keyword,comment,&status);
	else status = 0;
	if (fits_read_key(lodiffuseFits,TSTRING,"STAN_CON",keyword,comment,&status) == 0)
		fits_update_key(outFits,TSTRING,"STAN_CON",keyword,comment,&status);
	else status = 0;
	if (fits_read_key(lodiffuseFits,TSTRING,"FILE_ID",keyword,comment,&status) == 0)
		fits_update_key(outFits,TSTRING,"FILE_ID",keyword,comment,&status);
	else status = 0;

	/**
	long nrows = ExpMap.GetNrows();
	long ncols = ExpMap.GetNcols();
	*/
	long nrows = ExpMap.Rows();
	long ncols = ExpMap.Cols();

	long lodiffusepixel[3] = {1,1,1};
	long hidiffusepixel[4] = {1,1,1,1};

	bool inhiresmap = FALSE;
//	const double loconversion = fabs(lodl * D2R * lodb * D2R);
//	const double hiconversion = fabs(hidl * D2R * hidb * D2R);
	double ll, bb;
	
	for (outpixel[0]=1;outpixel[0]<=nrows;outpixel[0]++)
		for(outpixel[1]=1;outpixel[1]<=ncols;outpixel[1]++){
			ll = ExpMap.l(outpixel[0]-1,outpixel[1]-1);
			bb = ExpMap.b(outpixel[0]-1,outpixel[1]-1);
//			outtextstream << outpixel[0] << " " << outpixel[1] << " " << ll << " " << bb;
			if ( hiresstatus != 0 )  
				inhiresmap = FALSE;
			else {
//			    long rawpix = (long)(hix0 + (ll - hil0 + 0.5*hidl) / hidl)+1;
			    hidiffusepixel[0] = (long)(0.5 + fmod(hix0 + (ll - hil0) / hidl + 720.0/fabs(hidl), 360.0/fabs(hidl)));
			    hidiffusepixel[1] = (long)(0.5 + hiy0 + (bb - hib0) / hidb);
			    inhiresmap = hidiffusepixel[0] >= 1 && hidiffusepixel[0] <= hinx && hidiffusepixel[1] >= 1 && hidiffusepixel[1] <= hiny;
//			    outtextstream << " " << hidiffusepixel[0] << " " << hidiffusepixel[1] ;
			}
			if ( inhiresmap ) {
			    fits_read_pix(hidiffuseFits, TDOUBLE, hidiffusepixel, 1, &nulval, &pixels, &anynul, &hiresstatus) ;
				if ( hiresstatus == 0 && anynul == 0 ) {
					fits_write_pix(outFits, TDOUBLE, outpixel, 1, &pixels, &status) ; 
//					outtextstream << " " << pixels ; 
				}
				else {
//					fits_read_errmsg(errormessage);
//					fits_get_errstatus(hiresstatus, errorstatus);
//					outtextstream << " " << errormessage << " " << hiresstatus << " " << errorstatus ;
					hiresstatus = 0;
					inhiresmap = FALSE;
				}
			}
//			outtextstream << " " << inhiresmap ;
			if ( !inhiresmap ) {
//				long rawpix = (long)(lox0 + (ll - lol0 + 0.5*lodl) / lodl)+1;
				lodiffusepixel[0] = (long)(0.5 + fmod(lox0 + (ll - lol0) / lodl + 720.0/fabs(lodl), 360.0/fabs(lodl)));
				lodiffusepixel[1] = (long)(0.5 + loy0 + (bb - lob0) / lodb);
				if (lodiffusepixel[1] > lony)
					lodiffusepixel[1] = lony;
				else if (lodiffusepixel[1] < 1)
					lodiffusepixel[1] = 1;
//				outtextstream << " " << lodiffusepixel[0] << " " << lodiffusepixel[1] ;
				
				fits_read_pix(lodiffuseFits, TDOUBLE, lodiffusepixel, 1, NULL, &pixels, NULL, &status) ;
				fits_write_pix(outFits, TDOUBLE, outpixel, 1, &pixels, &status) ; 
//				outtextstream << " " << pixels ; 
			}
//			outtextstream << endl;
			if (status) throw;
		}
		
//	fits_write_2d_dbl(outFits, 0, ncols, nrows, ncols, outmap, &status);

	
	
// 	outtextstream.close();
	fits_close_file(outFits, &status);
	fits_close_file(expFits, &status);
	fits_close_file(lodiffuseFits, &status);
	fits_close_file(hidiffuseFits, &hiresstatus);

	return status;
	}


 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char expfile[FLEN_FILENAME];
	char loresdiffusefile[FLEN_FILENAME];	
	char hiresdiffusefile[FLEN_FILENAME] = "";	
	char outfile[FLEN_FILENAME];

	
	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("expfile", expfile);
	status = PILGetString("diffusefile", loresdiffusefile);	
	status = PILGetString("hiresdiffusefile", hiresdiffusefile);	
	status = PILGetString("outfile", outfile);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_gasmapgen5.cpp v.2 - 10/10/10 - A.C., A.T. ########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Enter exposure file name = "<<expfile<< endl;
	cout << "Enter output file name = " <<outfile<< endl;
	cout << "Enter diffuse model file name = "<< loresdiffusefile << endl;
	cout << "Enter high res diffuse model file name = "<< hiresdiffusefile << endl;
	cout << " "<< endl;
	cout << " "<< endl;	
	

	cout << "AG_gasmapgen5...............................starting"<< endl;		
	if (status == 0)	
		status = AG_gasmapgen(expfile, outfile, loresdiffusefile, hiresdiffusefile);
	cout << "AG_gasmapgen5............................... exiting"<< endl;		
	if (status) {
		if (status != 105) {
			if (outfile[0] == '!')
				remove(outfile+1);
			else
				remove(outfile);
			}
		printf("AG_gasmapgen5..................... exiting AG_gasmapgen ERROR:");		
		fits_report_error(stdout, status);	
		return status;			
		}			
	else {
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_gasmapgen........... exiting #################\n");
		printf("#################################################################\n\n\n");					
		}			
	
	return status;
}

