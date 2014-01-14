////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       Alike
//       Release: V0.0 -  26/Gen/2005
//       Contributors: Andrew Chen, Andrea Giuliani, Stefano Vercellone, Alberto Pellizzoni,
//		Alessio Trois, Sandro Mereghetti (IASF-Milano)
//	 V1.2 - 31 Jul 2009
//	 Contributors: Andrew Chen, Tomaso Contessi
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       26/Gen/2005
//                      First release: V1.0
//       		Author: Andrew Chen, Alessio Trois (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////



#include "RoiMulti5.h"
#include "PilParams.h"

using namespace std;

const PilDescription c_params[] = {
	{ PilString, "maplist", "Map file name" },
	{ PilString, "sarfile", "SAR file name" },
	{ PilString, "edpfile", "EDP file name" },
	{ PilString, "psdfile", "PSD file name" },
//	{ PilString, "expcorrfile", "Exposure Correction file name" },
	{ PilReal,   "ranal", "Radius of analysis" },
	{ PilInt,    "galmode", "Galactic parameter mode" },
	{ PilInt,    "isomode", "Isotropic parameter mode" },
	{ PilString, "srclist", "Sources list" },
	{ PilString, "extsrclist", "Extended sources list" },
	{ PilString, "outfile", "Output file name prefix" },
	{ PilReal,   "ulcl",    "Upper limit confidence level" },
	{ PilReal,   "loccl",   "Location contour confidence level" },
	{ PilNone,   "",   "" }
	};


class MultiExtParams: public PilParams
{
public:
	MultiExtParams(): PilParams(c_params) {}
	void Print() { PilParams::Print(c_params); }
};



/**

static int DoPIL(
	int argc,
	char *argv[],
	char* maplistname,
	char* sarfilename,
	char* edpfilename,
	char* psdfilename,
	char* expcorrfilename,
	char* srcfilename,
	char* extsrcfilename,
	char* outfilename,
	double& ranal,
	double& ulcl,
	double& loccl,
	int& galmode,
	int& isomode)
{
int status = PILInit(argc, argv);
int numpar = 0;
status = status || PILGetNumParameters(&numpar);
GetPilParam("maplist", maplistname, &status);

GetPilParam("sarfile", sarfilename, &status);
GetPilParam("edpfile", edpfilename, &status);
GetPilParam("psdfile", psdfilename, &status);
GetPilParam("expcorrfile", expcorrfilename, &status);

GetPilParam("ranal", &ranal, &status);
GetPilParam("galmode", &galmode, &status);
GetPilParam("isomode", &isomode, &status);
GetPilParam("srclist", srcfilename, &status);
GetPilParam("extsrclist", extsrcfilename, &status);
GetPilParam("outfile", outfilename, &status);
GetPilParam("ulcl", &ulcl, &status);
GetPilParam("loccl", &loccl, &status);

PILClose(status);
return status;
}


static void PrintDiffMode(int mode)
{
cout << mode;
if (mode==0)
	cout << " (constant)" << endl;
else if (mode==1)
	cout << " (from file)" << endl;
else if (mode==2)
	cout << " (variable)" << endl;
else if (mode==3)
	cout << " (single variable parameter)" << endl;
}



static void PrintInput(
	const char* maplistname,
	const char* sarfilename,
	const char* edpfilename,
	const char* psdfilename,
	const char* expcorrfilename,
	const char* srcfilename,
	const char* extsrcfilename,
	const char* outfilename,
	double ranal,
	double ulcl,
	double loccl,
	int galmode,
	int isomode)
{
cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;

cout << "Map file name = " << maplistname << endl;
cout << "SAR file name : "<< sarfilename << endl;
cout << "EDP file name : "<< edpfilename << endl;
cout << "PSD file name : "<< psdfilename << endl;
cout << "Exposure Correction file name : "<< expcorrfilename << endl;
cout << "Radius of analysis : " << ranal << endl;

cout << "Galactic parameter mode : "; PrintDiffMode(galmode);
cout << "Isotropic parameter mode : "; PrintDiffMode(isomode);

cout << "Sources list : " <<  srcfilename << endl;	
cout << "Extended sources list : " <<  extsrcfilename << endl;	
cout << "Output file name prefix : " <<  outfilename << endl;

cout << "Upper limit confidence level: " << ulcl << endl;
cout << "Location contour confidence level: " << loccl << endl;

cout << endl << endl;
}

*/

class AppScreen
{
public:
	AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "###### AG_multiExt  v.1.0 - 30/08/2011 - A.C., T.C., A.T. #########"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	}

	~AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "##########  Task AG_multiExt......... exiting ###################"<< endl;
	cout << "#################################################################"<< endl;
	}
};



int main(int argc, char *argv[])
{
AppScreen appScreen;

MultiExtParams mPars;
if (!mPars.Load(argc, argv))
	return -1;
cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;
mPars.Print();

const char* maplistname = mPars["maplist"];
const char* sarfilename = mPars["sarfile"];
const char* edpfilename = mPars["edpfile"];
const char* psdfilename = mPars["psdfile"];
// const char* expcorrfilename = mPars["expcorrfile"]; /// To Be Removed
double ranal = mPars["ranal"];
int galmode = mPars["galmode"];
int isomode = mPars["isomode"];
const char* srcfilename = mPars["srclist"];
const char* extsrcfilename = mPars["extsrclist"];
const char* outfilename = mPars["outfile"];
double ulcl = mPars["ulcl"];
double loccl = mPars["loccl"];

MapList maplist;
int mapCount = maplist.Read(maplistname);
if (!mapCount) {
	cerr << "File " << maplistname << " missing or empty" << endl;
	return -1;
	}

MapData mapData;
if (!mapData.Load(maplist))
	return -1;

SourceDataArray srcArr = ReadSourceFile(srcfilename);
if (!srcArr.Count())
	cout << "Warning: no sources loaded" << endl;

ExtList extList;
if (!extList.Read(extsrcfilename))
	return -1;

ExtData extData;
if (!extData.Load(extList))
	return -1;

RoiMulti roiMulti;
if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename))
	return -1;
if (!roiMulti.SetMaps(mapData , galmode, isomode))
	return -1;
if (!roiMulti.SetExtendedSources(extData))
	return -1;

string fileName;

fileName = outfilename;
fileName += ".log";
roiMulti.SetLogfile(fileName.c_str());
if (roiMulti.DoFit(srcArr, ranal, ulcl, loccl, 1))
	return -1;

roiMulti.Write(outfilename);
roiMulti.WriteSources(outfilename);
roiMulti.WriteHtml(outfilename);

fileName = outfilename;
fileName += "2";
roiMulti.Write(fileName.c_str(), false);

return 0;
}
