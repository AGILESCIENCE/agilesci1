
/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 *
 * Usage:
 *
 * ./bin/AG_add-ring-to-aitof-map
 *       <input aitof map>
 *       <input ring>
 *       <erase> -> if yes/y the pixels values of 'input aitof map' will be set to 0 at the beginning.
 *       <destName [OPTIONAL]>  -> if not set, the 'the input aitof map' will be overwritten. If set another file will be created.
 *
*/


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <cstdio>

#include <PilParams.h>
#include <prj.h>
#include "AgileMap.h"
#include "FitsUtils.h"


using std::cout;
using std::endl;
using std::cerr;

enum ProjectionType {ARC, AIT};

const char* startString = {
"################################################################\n"
"###             Task AG_add-ring-to-aitof-map v0.1           ###"
};

const char* endString = {
"### Task AG_add-ring-to-aitof-map exiting ..................###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
	{ PilString,"mapPath","Map file path"},
 	{ PilString,"ringPath","Ring file path"},
 	{ PilBool,"eraseOldData","Erase data in input map? [yes/y, np/n]"},
	{ PilString,"destinationName","Output map name. If null the input map will be overwritten"},
	{ PilNone,"",""}
};

void printMapInfo(AgileMap& inputMap){
	cout << "\nFilename: " << inputMap.GetFileName() << endl;
	cout << "InputMap dims: " << inputMap.Dims() << endl;
	cout << "InputMap rows: " << inputMap.Rows() << endl;
	cout << "InputMap dim(0): " << inputMap.Dim(0) << endl;
	cout << "InputMap cols: " << inputMap.Cols() << endl;
	cout << "InputMap dim(1): " << inputMap.Dim(1) << endl;
  cout << "InputMap size: " << inputMap.Size() << endl;
	cout << "(L,B): " << inputMap.GetMapCenterL() << " , " << inputMap.GetMapCenterB() << endl;
	cout << "(X0,Y0): " << inputMap.GetX0() << " , " << inputMap.GetY0() << endl;
	cout << "(Xbin,Ybin): " << inputMap.GetXbin() << " , " << inputMap.GetYbin() << endl;
}

inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

bool removeIfExist(std::string filename){
	if( exists(filename) )
	{
		cout << "The file " << filename << " already exists..removing it.." << endl;
		remove(filename.c_str());
		return true;
	}
	return false;
}
bool inmap(int i, int ii, int mxdim)
{
    return (i < mxdim) && (i >= 0) &&
           (ii < mxdim) && (ii >= 0);
}


/*
AITOFF:
CDELT1  = -0.5
CDELT2  =  0.5
CUNIT1  = 'deg     '
CUNIT2  = 'deg     '

RING:
CDELT1  =  -0.1
 CDELT2  =  0.1

 la=ba=0
 mdim = 360 = CRPIX1
*/

int main(int argc, char*argv[]){

	cout << startString << "\n" << endl;


	/////////////////////////////////////////////////////////////////
	// INPUT PARAMETERS
	PilParams params(paramsDescr);

  if (!params.Load(argc, argv))
      return EXIT_FAILURE;

	const char * mapPath 	  = params["mapPath"];
	const char * ringPath   = params["ringPath"];
	bool eraseOldData			  = params["eraseOldData"];
	const char * destinationName = params["destinationName"];

	std::string _mapPath(mapPath);
	std::string _ringPath(ringPath);
	std::string _destinationName(destinationName);

	params.Print();

	if(!exists(mapPath))
	{
		cerr << "\nThe input file " << mapPath << " does not exist!" <<endl;
		exit (EXIT_FAILURE);
	}
	if(!exists(ringPath))
	{
		cerr << "\nThe input file " << ringPath << " does not exist!" <<endl;
		exit (EXIT_FAILURE);
	}

	AgileMap inputMap(mapPath);
	printMapInfo(inputMap);

	AgileMap ring(ringPath);
	printMapInfo(ring);

	cout << endl;

	if(eraseOldData)
	{
		cout << "Erasing data of " << mapPath << ".." << endl;
		inputMap.Zero();
	}

	/////////////////////////////////////////////////////////////////
	// PROJECTION PARAMETERS

	struct prjprm *prj = new(prjprm);
	prjini(prj);
	ProjectionType proj = AIT;
	double laa = inputMap.GetMapCenterL() * DEG2RAD;
	double baa = inputMap.GetMapCenterB() * DEG2RAD;
	double mres = fabs(inputMap.GetXbin());							/// CDELT1
	double mdim = inputMap.GetX0();                     /// CRPIX1
	long mxdim = inputMap.Rows(); // dimension (in pixels) of the map


	/////////////////////////////////////////////////////////////////
	// CORE ALGORITHM

	double z = 0, zz = 0;
	for(int i=0; i<ring.Rows(); i++){
		for(int j=0; j<ring.Cols(); j++){

			int ringBin = ring.Bin(i,j);
			double ringValue = ring[ringBin];
			double ringL = ring.l(i,j);
			double ringB = ring.b(i,j);

			double l = ringL*DEG2RAD;
			double b = ringB*DEG2RAD;

			double theta = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (theta < -1.0)
					theta = M_PI;
			else if (theta > 1.0)
					theta = 0.0;
			else
					theta = acos(theta);

			l=l-laa;

			if (l < M_PI)
				l=-l;
			else
				l=2*M_PI -l;

 			double mapColAit=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) );
			double mapRowAit=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

			z=(int)floor(((mapColAit+(mdim/2.))/mres));
			zz=(int)floor(((mapRowAit+(mdim/2.))/mres));

  		int mapBin = inputMap.Bin(z,zz);
			double mapValue = inputMap[mapBin];

			if (inmap(z, zz, mxdim) && mapValue < ringValue) {
					inputMap[mapBin] = ringValue;

			}
		}
	}


	/////////////////////////////////////////////////////////////////
	// I/O
	cout << endl;

	std::string _tempDestinationName("abcde.gz");

	FitsFile f;
	if (!f.Create(_tempDestinationName.c_str())) {
		cerr << "ERROR " << f.Status() << " creating " << _tempDestinationName << endl;
		exit (EXIT_FAILURE);
	}

	MatD mat;
	inputMap.TransposeTo(mat);

 	FitsFile from;
	if (!from.Open(mapPath)) {
	    cerr << "ERROR " << f.Status() << " reopening " << inputMap.GetFileName() << endl;
	    exit (EXIT_FAILURE);
	}

	int status = 0;
	fits_copy_header(from, f, &status);
	if (status) {
	    cerr << "ERROR " << status << " copying hdu from " << inputMap.GetFileName() << " to " << _tempDestinationName << endl;
	    exit (EXIT_FAILURE);
	}
	from.Close();

	f.UpdateKey("NAXIS", 2);
	f.UpdateKey("NAXIS1", mat.Cols());
	f.UpdateKey("NAXIS2", mat.Rows());

	long fpixel[2] = { 1, 1 };
	fits_write_pix(f, TDOUBLE, fpixel, mat.Size(), const_cast<double*>(mat.Buffer()), f);

	f.UpdateKey("CTYPE1", "GLON-AIT");
	f.UpdateKey("CTYPE2", "GLAT-AIT");
	f.UpdateKey("CRPIX1", inputMap.GetX0());
	f.UpdateKey("CRVAL1", inputMap.GetMapCenterL());
	f.UpdateKey("CDELT1", inputMap.GetXbin());
	f.UpdateKey("CRPIX2", inputMap.GetY0());
	f.UpdateKey("CRVAL2", inputMap.GetMapCenterB());
	f.UpdateKey("CDELT2", inputMap.GetYbin());
	f.UpdateKey("LONPOLE", inputMap.GetLonpole());
	f.UpdateKey("MINENG", inputMap.GetEmin());
	f.UpdateKey("MAXENG", inputMap.GetEmax());
	f.UpdateKey("INDEX", inputMap.GetMapIndex());
	f.UpdateKey("SC-Z-LII", inputMap.GetLpoint());
	f.UpdateKey("SC-Z-BII", inputMap.GetBpoint());
	// f.UpdateKey("SC-LONPL", inputMap. m_gp); --> no getter is available

	if (inputMap.GetStartDate()[0])
		f.UpdateKey("DATE-OBS", inputMap.GetStartDate());
	if (inputMap.GetEndDate()[0])
		f.UpdateKey("DATE-END", inputMap.GetEndDate());

	f.UpdateKey("TSTART", inputMap.GetTstart());
	f.UpdateKey("TSTOP", inputMap.GetTstop());
	f.UpdateKey("FOVMIN", inputMap.GetFovMin());
	f.UpdateKey("FOV", inputMap.GetFovMax());
	f.UpdateKey("ALBEDO", inputMap.GetAlbedo());
	f.UpdateKey("PHASECOD", inputMap.GetPhaseCode());

	if (inputMap.GetStep())
		f.UpdateKey("STEP", inputMap.GetStep());

	if (inputMap.GetSkyL()[0])
		f.UpdateKey("SKYL", inputMap.GetSkyL());
	if (inputMap.GetSkyH()[0])
		f.UpdateKey("SKYH", inputMap.GetSkyH());

	if (f.Status())
	{
		cerr << "ERROR " << f.Status() << " writing to " << _tempDestinationName << endl;
		exit (EXIT_FAILURE);
	}
	f.Close();


	// overwrite inputMap?
	if(_destinationName.compare("null") == 0 )
	{
		_destinationName = _mapPath;
	}

	if( exists(_destinationName) )
	{
		remove(_destinationName.c_str());
	}
	cout << "Creating " << _destinationName << " fits file.. " << endl;
	rename(_tempDestinationName.c_str(),_destinationName.c_str());


	cout << endl << endString << endl;
	return 0;
}
