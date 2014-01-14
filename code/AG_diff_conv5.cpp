////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_diff_conv
//       Release: V2.0 -  17/Oct/2010
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
//       17/Oct/2010
//                      Second release: V2.0
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;



#include <cstdio>
#include <cstring>
#include <cmath>
//#include <fstream>

#include "wcsmath.h"
#include "wcstrig.h"
#include "fitsio.h"
#include "pil.h"

#include "CalibUtils.h"
#include "MathUtils.h"


class AlikeNoDispPsfArray: public AeffGrid, public PsfGrid, public EdpGrid
{
public:
	AlikeNoDispPsfArray(const char* psfName, const char* raeffName, const char* edpName):
		AeffGrid(raeffName), PsfGrid(psfName), EdpGrid(edpName) {}
	~AlikeNoDispPsfArray() {}
public:
	VecD CalcWeighted(float index, float E_low, float E_high, float theta);
};



VecD AlikeNoDispPsfArray::CalcWeighted(float index, float E_low, float E_high, float theta)
{
if (index < 0)
	index = -index;

const VecF& psfenergies = PsfGrid::Energies();
int numpsfenergies = psfenergies.Size();
const VecF& rhos = PsfGrid::Rhos();
int numrhos = rhos.Size();
double psfinc = rhos[1]-rhos[0];
cout << "Num PSF energies = " << numpsfenergies << endl;
MatD psftable(numrhos, numpsfenergies);
VecD ringsr(numrhos);

int i=0, j=0, etrue=0, eobs=0;
for (i=0; i<numrhos; ++i)
	ringsr[i] = 2.0 * PI * psfinc * sind(rhos[i]);

for (j=0; j<numpsfenergies; ++j) {
	double totaltemp = 0.0;
	for (i=0; i<numrhos; ++i) {
		psftable(i, j) = PsfGrid::Val(rhos[i], 0, theta, 0, psfenergies[j]) / sind(rhos[i]);
//		cout << i << " " << j << " " << rhos[i] << " " << psfenergies[j] << " " << psftable(i,j)*ringsr[i]/sind(rhos[i]) << endl;
		totaltemp += psftable(i,j) * ringsr[i] / sind(rhos[i]);
		}
	for (i=0; i<numrhos; ++i)
		psftable(i, j) /= totaltemp;
	}

int lowetrue = psfenergies.GeomIndex(E_low);
int highetrue = psfenergies.GeomIndex(E_high);
//cout << "Low E index = " << lowetrue << ", High E index = " << highetrue << endl;


MatD specwt(numpsfenergies,1);
MatD edparr(numpsfenergies,1);
specwt = 0.0;
edparr = 0.0;

for (etrue = lowetrue; etrue <= highetrue; etrue++) {
	if (etrue < numpsfenergies - 1)
		specwt[etrue] = pow(psfenergies[etrue], 1.0f-index) - pow(psfenergies[etrue+1], 1.0f -index);
	else
		specwt[etrue] = pow(psfenergies[etrue], 1.0f-index);
	for (eobs = lowetrue;  eobs <= highetrue; eobs++)
		edparr[etrue] = edparr[etrue] + EdpGrid::Val(psfenergies[etrue], psfenergies[eobs], theta, 0.0f);
//	cout << psfenergies[etrue] << " " << specwt[etrue] << " " << edparr[etrue] << endl;
	}
	
MatD monoArrays(numrhos, numpsfenergies);
for (etrue = lowetrue; etrue <= highetrue; etrue++)
	for (i = 0; i < numrhos; i++) {
		monoArrays(i, etrue) = 
				PsfGrid::Val(rhos[i], 0.0f, theta, 0.0f, psfenergies[etrue]) *
				AeffGrid::Val(psfenergies[etrue], theta, 0.0f) *
				edparr[etrue];
//		cout << psfenergies[etrue] << " " << rhos[i] << " " << monoArrays(i, etrue) << endl;
	}


MatD psfMat = monoArrays * specwt;
VecD psfArray(psfMat.Size());
	psfArray = psfMat.Buffer();

double deltaE2 = (highetrue-lowetrue) * (highetrue-lowetrue);
psfArray /= deltaE2;
double dsum = 0.0;
for (i=0; i<numrhos; i++)
	dsum += psfArray[i] * ringsr[i];
psfArray /=  dsum;
return psfArray;
}








double psflookup(int m_psfCount, const VecF& m_rhoArr, const VecD& m_psfArr, double srcdist)
{
double deltaRho = m_rhoArr[1]-m_rhoArr[0];
double maxRho = m_rhoArr[m_psfCount-1]+deltaRho/2;
if (srcdist>=maxRho)
	return 0;
else {
	int psfind = (int)(srcdist/deltaRho);
	double resid = (srcdist - m_rhoArr[psfind])/deltaRho;
	double correction;
	if (resid<0)
		if (psfind==0)
			correction = m_psfArr[1]-m_psfArr[0];
		else
			correction = m_psfArr[psfind]-m_psfArr[psfind-1];
		else
			if (psfind==m_psfCount-1)
				correction = -m_psfArr[psfind];	/// Next one would be zero
			else
				correction = m_psfArr[psfind+1]-m_psfArr[psfind];
	return m_psfArr[psfind] + resid * correction;
	}
}

int AG_diffuse_convolve(char * diffusefile, char * psdfile, char * sarfile, char * edpfile, char *  outfile){
	
	int status = 0;

	char logfile[FLEN_FILENAME];
	strcpy(logfile,outfile);
	strcat(logfile,".log");
//	ofstream ofout(logfile);
	fitsfile * diffuseFits;
	if ( fits_open_file(&diffuseFits, diffusefile, READONLY, &status) != 0 ) {
		printf("Errore in apertura file '%s'\n", diffusefile);
		return status;
	}
	else {
		
		fits_movabs_hdu(diffuseFits, 2, NULL, &status);	
		
		int bitpix   =  DOUBLE_IMG;
		int naxis;   
		long naxes[3] ;
		
		fits_get_img_param(diffuseFits, 3, &bitpix, &naxis, naxes, &status);
		cout << naxis << " axes: (" << naxes[0] << ", " << naxes[1] << ", " << naxes[2] << ")" << endl;
		
		double x0, y0, emin, emax;	
		double l0, b0;
		double dl, db;	
		double index;
		fits_read_key(diffuseFits,TDOUBLE,"CRVAL1",&l0,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"CDELT1",&dl,NULL,&status);
		double fdl = fabs(dl);
		fits_read_key(diffuseFits,TDOUBLE,"CRPIX1",&x0,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"CRVAL2",&b0,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"CDELT2",&db,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"CRPIX2",&y0,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"E_MIN",&emin,NULL,&status);
		fits_read_key(diffuseFits,TDOUBLE,"E_MAX",&emax,NULL,&status);
		int oldstatus = status;
		fits_read_key(diffuseFits,TDOUBLE,"INDEX",&index,NULL,&status);
		if (status != 0) {
			cout << "Spectral index not found; set to -2.0" << endl ;
			status = oldstatus;
			index = -2.0;
		}
		else if (index > 0)
			index = -index;
		
		
		double polecosb = 1.0 - sind(90.0 - db);
		long pixel[2] = { 1, 1 };
		double * diffuse2 = new double[naxes[0] * naxes[1]];
		fits_read_pix(diffuseFits, TDOUBLE, pixel, naxes[0]*naxes[1], NULL, diffuse2, NULL, &status);
		cout << "Diffuse array successfully created" << endl;
		
		
		cout << "l = " << l0 << " + " << dl << " * (x - " << x0 << ")" << endl;
		cout << "b = " << b0 << " + " << db << " * (y - " << y0 << ")" << endl;
		cout << "Emin = " << emin << ", Emax = " << emax << endl;
		
		AlikeNoDispPsfArray psfarray(psdfile, sarfile, edpfile);
		const VecF& rhoArr = psfarray.Rhos();
		VecD psfArr = psfarray.CalcWeighted(index, emin, emax, 30.0);
		int nrhos = rhoArr.Size();
//		for (int i=0 ; i< nrhos ; i++) {ofout << rhoArr[i] << " " << psfArr[i] << endl;}
		
		int height = 2 * (int)ceil(rhoArr[nrhos-1] / db + 0.5);
		
		cout << "Last angle = " << rhoArr[nrhos-1] << ", height = " << height << endl;
		
		fitsfile * outFits;
		if ( fits_create_file(&outFits, outfile, &status) != 0 ) {
			printf("Errore in apertura file '%s'\n", outfile);
			return status;
		}
		else {
			fits_copy_file(diffuseFits, outFits, 1, 0, 0, &status);
			
			double distance;
			long ylow, yhigh;
			
			long jlow = 0;
			long jhigh = naxes[1];
			cout << "Old array indices from " << jlow << " to " << jhigh << endl;
			cout << "Low b = " << b0+db*(1-y0) << " and high b = " << b0+db*(jhigh-y0) << endl;
			double nearpole = 90 - 4 * fabs(db);
			if (b0+db*(1-y0) > -nearpole)
				jlow += height / 2;
			if (b0+db*(jhigh-y0) < nearpole)
				jhigh -= height / 2;
			cout << "New array indices from " << jlow << " to " << jhigh << endl;
			long ilow = 0;
			long ihigh = naxes[0];
			if (fabs(naxes[0] - 360.0 / fdl) > 4) {
				int halfwidth = (2 * (int)ceil(rhoArr[nrhos-1] / fdl + 0.5)) / 2;
				ilow = halfwidth;
				ihigh = naxes[0] - halfwidth;
			}
			long maskwidth = long(0.5 + 360.0 / fdl);
			double * mask2 = new double[maskwidth * height];
			
			long numnewgas2 = (ihigh-ilow) * (jhigh-jlow);
			double * newgas2 = new double[numnewgas2];
			for (long j = 0; j < numnewgas2 ; j++)
				newgas2[j] = 0;
			
			double l00 = l0+dl*(1-x0);
			for (long j = 0; j < naxes[1] ; j++) {
				double b00 = b0+db*(j+1-y0);
				if (b00 >= 90.0) b00 = 90.0;
				if (b00 <= -90.0) b00 = -90.0;
				ylow = j - height / 2;
				if (ylow < 0) {ylow = 0;}
				yhigh = j + height / 2;
				if (yhigh > naxes[1]) {yhigh = naxes[1];}
				double masktot = 0;
				cout << "j = " << j << ", (l00, b00) = (" << l00 << "," << b00 << ")" << endl;
				cout << "ylow = " << ylow << ", yhigh = " << yhigh << ", blow = " << b0+db*(ylow+1-y0) << ", bhigh = " << b0+db*(yhigh+1-y0) << endl; 
				for (long m = ylow; m <  yhigh ; m++) {
					double b = b0+db*(m+1-y0);
					double cosb = 0.0;
					if (b >= 90.0) {
						b = 90.0;
						cosb = polecosb;
					} else if (b <= -90.0) {
						b = -90.0;
						cosb = polecosb;
					} else 
						cosb = fabs(sind(b+0.5*db)-sind(b-0.5*db));
					for (long k = 0; k < maskwidth ; k++) {
						distance = SphDistDeg(l00, b00, l0+dl*(k+1-x0), b);
//						ofout << l00 << " " << b00 << " " << l0+dl*(k+1-x0) << " " << b << " " << distance << endl;
						if (distance > height / 2) {
							mask2[k + (m-ylow) * maskwidth] = 0;
						} else {
							mask2[k + (m-ylow) * maskwidth] = psflookup(nrhos, rhoArr, psfArr, distance) * cosb;
//							ofout << k << " " << m << " " << distance << " " << mask2[k + (m-ylow) * maskwidth] << endl; 
							masktot += mask2[k + (m-ylow) * maskwidth];
						}
					}			
				}
				for (long m = ylow; m <  yhigh ; m++) {
					for (long k = 0; k < maskwidth ; k++) {
						mask2[k + (m-ylow) * maskwidth] *= 1.0/masktot;
//						ofout << m - ylow << " " << k << " " << mask2[k + (m-ylow) * maskwidth] << endl;
					}
				}

/*				fitsfile *maskfits;
				char maskname[FLEN_FILENAME];
				char dummy[FLEN_FILENAME];
//				strcpy(maskname, "!");
				strcpy(maskname, outfile);
				strcat(maskname, "_");
				sprintf(dummy, "%li", j);
				strcat(maskname, dummy);
				strcat(maskname, ".gz");
				fits_create_file(&maskfits, maskname, &status);
				cout << "creating file " << maskname << ", status = " << status << endl;
				long maskaxes[2] = {maskwidth, yhigh-ylow};
				long maskfpixel[2] = {1,1};
				fits_create_img(maskfits, DOUBLE_IMG, 2, maskaxes, &status);
				cout << "creating image " << maskname << ", status = " << status << endl;
				fits_write_pix(maskfits, TDOUBLE, maskfpixel, maskwidth * (yhigh-ylow), mask2, &status);
				cout << "writing image " << maskname << ", status = " << status << endl;
				fits_close_file(maskfits, &status);
				cout << "closing file " << maskname << ", status = " << status << endl;
*/
				long mlow = (jlow < ylow) ? ylow : jlow;
				long mhigh = (jhigh < yhigh) ? jhigh : yhigh;
				cout << "mlow = " << mlow << ", mhigh = " << mhigh << endl;
				for (long i = 0; i < naxes[0] ; i++) {
					double pixval = diffuse2[i + j * naxes[0]]; 
//					ofout << i << " " << pixval << endl;
					for (long m = mlow; m < mhigh ; m++) {
						for (long k = ilow ; k < ihigh ; k++) {
							newgas2[(k - ilow) + (m-jlow) * (ihigh - ilow)] += pixval * mask2[ (k > i ? k - i : i - k) + (m - ylow) * maskwidth];
//							ofout << m << " " << k << " " << (k > i ? k - i : i - k) << " " << m - ylow << " " << mask2[ (k > i ? k - i : i - k) + (m - ylow) * maskwidth] << endl;
						}
						/* for (long k = i; k < naxes[0] ; k++) {
							if (k >= ilow && k < ihigh)
								newgas2[(k - ilow) + (m-jlow) * (ihigh - ilow)] += pixval * mask2[ (k - i) + (m - ylow) * maskwidth];
						}
						for (long k = 0; k < i ; k++) {
							if (k >= ilow && k < ihigh)
								newgas2[(k - ilow) + (m-jlow) * (ihigh - ilow)] += pixval * mask2[ (i - k) + (m - ylow) * maskwidth];
						}
						 */
					}
				}
			}
			long outnaxes[2] = {ihigh-ilow,jhigh-jlow};
			char keywordstring[FLEN_KEYWORD];
			char comment[FLEN_COMMENT];
			
			fits_create_img(outFits, FLOAT_IMG, 2, outnaxes, &status);
			fits_write_pix(outFits, TDOUBLE, pixel, (ihigh-ilow)*(jhigh-jlow), newgas2, &status) ;
			
//			l0 += dl * ilow;
			fits_update_key(outFits,TDOUBLE,"CRVAL1",&l0,NULL,&status);
			x0 -= ilow;
			fits_update_key(outFits,TDOUBLE,"CRPIX1",&x0,NULL,&status);
			fits_update_key(outFits,TDOUBLE,"CDELT1",&dl,NULL,&status);
			fits_read_key(diffuseFits,TSTRING,"CTYPE1",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"CTYPE1",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(diffuseFits,TSTRING,"CUNIT1",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"CUNIT1",keywordstring,comment,&status);
			else
				status = 0;
//			b0 += db * jlow;
			fits_update_key(outFits,TDOUBLE,"CRVAL2",&b0,NULL,&status);
			y0 -= jlow;
			fits_update_key(outFits,TDOUBLE,"CRPIX2",&y0,NULL,&status);
			fits_update_key(outFits,TDOUBLE,"CDELT2",&db,NULL,&status);
			fits_read_key(diffuseFits,TSTRING,"CTYPE2",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"CTYPE2",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(diffuseFits,TSTRING,"CUNIT2",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"CUNIT2",keywordstring,comment,&status);
			else
				status = 0;
			int pixcent = TRUE;
			fits_read_key(diffuseFits,TLOGICAL,"PIXCENT",&pixcent,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TLOGICAL,"PIXCENT",&pixcent,comment,&status);
			else
				status = 0;
			fits_read_key(diffuseFits,TSTRING,"BUNIT",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"BUNIT",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(diffuseFits,TSTRING,"PRIMTYPE",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"PRIMTYPE",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(diffuseFits,TSTRING,"EXTNAME",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"EXTNAME",keywordstring,comment,&status);
			else
				status = 0;
			fits_update_key(outFits,TDOUBLE,"E_MIN",&emin,"Energy range",&status);
			fits_update_key(outFits,TDOUBLE,"E_MAX",&emax,"Energy range",&status);
			fits_update_key(outFits,TDOUBLE,"INDEX",&index,NULL,&status);
			
			fitsfile * psdFits;
			fits_open_file(&psdFits, psdfile, READONLY, &status);
			fits_read_key(psdFits,TSTRING,"TELESCOP",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"TELESCOP",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(psdFits,TSTRING,"INSTRUME",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"INSTRUME",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(psdFits,TSTRING,"DH_CONF_",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"DH_CONF_",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(psdFits,TSTRING,"STAN_CON",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"STAN_CON",keywordstring,comment,&status);
			else
				status = 0;
			fits_read_key(psdFits,TSTRING,"FILE_ID",keywordstring,comment,&status);
			if (status == 0)
				fits_update_key(outFits,TSTRING,"FILE_ID",keywordstring,comment,&status);
			else
				status = 0;
			fits_close_file(psdFits, &status);
			fits_close_file(outFits, &status);
delete [] mask2;
			delete [] newgas2;
		}
		
		// clean up
		fits_close_file(diffuseFits, &status);
		delete [] diffuse2;
		
	}
//	ofout.close();
	return status;
}

 	
int main(int argc,char **argv)
{

	int status = 0, numpar = 0;
	char psdfile[FLEN_FILENAME];
	char sarfile[FLEN_FILENAME];
	char edpfile[FLEN_FILENAME];
	char diffusefile[FLEN_FILENAME];	
	char outfile[FLEN_FILENAME];

	status = PILInit(argc,argv);
	status = PILGetNumParameters(&numpar);
	status = PILGetString("diffusefile", diffusefile);	
	status = PILGetString("psdfile", psdfile);
	status = PILGetString("sarfile", sarfile);
	status = PILGetString("edpfile", edpfile);
	status = PILGetString("outfile", outfile);

	status = PILClose(status);

	cout << " "<< endl;
	cout << " "<< endl;	
	cout << "#################################################################"<< endl;
	cout << "########## AG_diffuse_convolve.cpp v.0 - 19/12/05 - A.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << " "<< endl;
	cout << "INPUT PARAMETERS:"<< endl;
	cout << " "<< endl;
	cout << "Diffuse model file name = "<< diffusefile << endl;
	cout << "Point spread function file name = "<<psdfile<< endl;
	cout << "Sensitive area file name = "<<sarfile<< endl;
	cout << "Energy dispersion file name = "<<edpfile<< endl;
	cout << "Enter output file name = " <<outfile<< endl;
	cout << " "<< endl;
	cout << " "<< endl;	
	

	cout << "AG_diffuse_convolve...............................starting"<< endl;		
	if (status == 0)	
		status = AG_diffuse_convolve(diffusefile, psdfile, sarfile, edpfile, outfile);
	cout << "AG_diffuse_convolve............................... exiting"<< endl;		
	if (status) {
/*		if (status != 105) {
			char * temp = new char[128];
			sprintf(temp, "rm %s", outfile);
			system(temp);
			}*/
		printf("AG_diffuse_convolve..................... exiting AG_diffuse_convolve ERROR:");		
		fits_report_error(stdout, status);	
		return status;			
		}			
	else {
		printf("\n\n\n###################################################################\n");
		printf("#########  Task AG_diffuse_convolve........... exiting #################\n");
		printf("#################################################################\n\n\n");					
		}			
	
	return status;
}

