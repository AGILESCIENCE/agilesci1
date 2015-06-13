


#include <cstring>
#include <fitsio.h>
#include <iostream>
#include <cstdlib>
#include "SkyMap.h"
#include "MathUtils.h"

using std::cerr;
using std::cout;
using std::endl;
//using std::isnan;

static bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, float defaultVal, int* status)
{
	if (!*status) {
		int status2 = 0;
		fits_read_key(fptr, TFLOAT, const_cast<char*>(name), value, NULL, &status2);
		if (status2==0)
			return true;
		else if (status2==KEY_NO_EXIST)
			*value = defaultVal;
		else
			*status = status2;
	}
	return false;
}

static void ReadFitsKey(fitsfile* fptr, const char* name, float* value, int* status)
{
	fits_read_key(fptr, TFLOAT, const_cast<char*>(name), value, NULL, status);
}


/* static void ReadFitsKey(fitsfile* fptr, const char* name, char* value, int* status)
{
	fits_read_key(fptr, TSTRING, const_cast<char*>(name), value, NULL, status);
}
 */
/// #define MORE_KEYWORDS

namespace
{
const char* const knownKeywords[] = {
	"XTENSION",
	"BITPIX",
	"NAXIS",
	"NAXIS1",
	"NAXIS2",
	"NAXIS3",
	"NAXIS4",
	"PCOUNT",
	"GCOUNT",
	"BSCALE",
	"BZERO",
	"CDELT1",
	"CRPIX1",
	"CRVAL1",
	"CDELT2",
	"CRPIX2",
	"CRVAL2",
	"EXTNAME",
#ifdef MORE_KEYWORDS
	"CTYPE1",
	"CUNIT1",
	"CTYPE2",
	"CUNIT2",
	"PIXCENT",
	"BUNIT",
	"PRIMTYPE",
	"INSTRUME",
#endif
	0};

bool IsUnknownKeyword(const char* keyword)
{
	for (int i=0; knownKeywords[i]; ++i)
		if (!strcmp(keyword, knownKeywords[i]))
			return false;
	return true;
}

void ReadUnknKeywords(fitsfile *fptr, int& count, int*& offsets, char*& text, int* status)
{
	count = 0;
	offsets = 0;
	text = 0;
	if (*status)
		return;
	
	int l_count;
	fits_get_hdrspace(fptr, &l_count, 0, status);
	if (l_count==0 || *status)
		return;
	char* l_text = new char[l_count*84];
	int* l_offsets = new int[l_count];
	int index = 0;
	int offset = 0;
	for (int i=0; i<l_count && !*status; ++i) {
		char keyname[128];
		char value[128];
		fits_read_keyn(fptr, i+1, keyname, value, 0, status);
		if (IsUnknownKeyword(keyname)) {
			char record[128];
			fits_read_record(fptr, i+1, record, status);
			int length = strlen(record);
			strcpy(l_text+offset, record);
			l_offsets[index] = offset;
			offset += length+1;
			++index;
		}
	}
	
	if (index && !*status) {
		count = index;
		offsets = new int[count];
		memcpy(offsets, l_offsets, sizeof(int)*count);
		text = new char[offset];
		memcpy(text, l_text, offset);
	}
	delete[] l_offsets;
	delete[] l_text;
}

}

SkyMap converttoSkyMap(const char* fileName)
	{
		SkyMap outmap;

		fitsfile *fptr;
		int status = 0;
		fits_open_image(&fptr, fileName, READONLY, &status);
		if (status)
			cerr << "Convert to SkyMap: failed opening file " << fileName << " Fits error: " << status << endl;
		else {
				int bitpix, naxis;
				long naxes[4];
				fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
				if (status)
					cerr << "Convert to SkyMap: Could not read the image parameters" << endl;
				if (status==0) {
					int nx = naxes[0];
					int ny = naxes[1];
					float* values = new float[nx*ny];
					for (int j=0; j<ny ; j++){
						long fpixel[4] = { 1, j+1, 1, 1};
						fits_read_pix(fptr, TFLOAT, fpixel, nx, 0, values+j*nx, 0, &status);
					}					
					if (status)
						cerr << "Convert to SkyMap: Error reading the image" <<  endl;
					float xRes, yRes, xRef, yRef, xPoint, yPoint;
					ReadFitsKey(fptr, "CDELT1", &xRes, &status);
					ReadFitsKey(fptr, "CDELT2", &yRes, &status);
					ReadFitsKey(fptr, "CRPIX1", &xRef, &status);
					ReadFitsKey(fptr, "CRPIX2", &yRef, &status);
					ReadFitsKey(fptr, "CRVAL1", &xPoint, &status);
					ReadFitsKey(fptr, "CRVAL2", &yPoint, &status);
					
					if (status)
						cerr << "Convert to SkyMap: Error reading the basic keywords" << endl;
					/* else {
						WarnCTYPE(fptr, fileName, "CTYPE1");
						WarnCTYPE(fptr, fileName, "CTYPE2");
						WarnCUNIT(fptr, fileName, "CUNIT1");
						WarnCUNIT(fptr, fileName, "CUNIT2");
						WarnCROTA(fptr, fileName, "CROTA1");
						WarnCROTA(fptr, fileName, "CROTA2");
					}
					 */					
					int kwrdCount;
					int* kwrdIndices;
					char* kwrdText;
					ReadUnknKeywords(fptr, kwrdCount, kwrdIndices, kwrdText, &status);
					/// PrintKeywords(fptr, &status);
					if (status)
						cerr << "Convert to SkyMap: Error reading ancillary keywords" << endl;
					else {

						int y_padwidth = 0;
						int new_ny = ny;
						float new_yRef = yRef; 
//						float b_lo = yPoint + yRes * (1 - new_yRef);
//						float b_hi = yPoint + yRes * (new_ny - new_yRef);
						int x_padwidth = 0;
						int new_nx = nx;
						float new_xRef = xRef;

						SkyMap map(new_ny, new_nx, xRes, yRes, xPoint, yPoint, new_xRef, new_yRef);
						for (int i=0; i<new_nx ; i++)
							for (int j=0; j<new_ny ; j++) {
								map.SetValue(j,i,0);
							}
						for (int i=0; i<nx ; i++)
							for (int j=0; j<ny ; j++) {
								if (isnan(values[i+j*nx]))
									values[i+j*nx]=0;
								map.SetValue(j+y_padwidth,i+x_padwidth,values[i+j*nx]);
						}
						
						char** keywords = new char*[kwrdCount+1];
						for (int k=0; k<kwrdCount; k++) {
							keywords[k] = new char[FLEN_CARD];
							strcpy(keywords[k],kwrdText+kwrdIndices[k]);
//							printf("%s\n",keywords[k]);
						}
						keywords[kwrdCount]=0;
						map.SetKeywords(keywords);
						for (int k=0; k<kwrdCount; k++) {
							delete [] keywords[k];
						}
						delete [] keywords;
						outmap = map;
				}
			}
			int closeStatus = 0;	/// Close this file in any case
			fits_close_file(fptr, &closeStatus);
		}
	
	return outmap;
}


int main(int argC, char* argV[])
{
if (argC<2) {
	cout << "Usage: AG_converttoSkyMap5 inputMap resultMapName" << endl;
	return 0;
	}


SkyMap map2 = converttoSkyMap(argV[1]);
bool saved = map2.Save(argV[2]);
return saved ? 0 : -1;
}

