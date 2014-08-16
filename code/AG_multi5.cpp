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
	{ PilString, "outfile", "Output file name prefix" },
	{ PilReal,   "ulcl",    "Upper limit confidence level" },
	{ PilReal,   "loccl",   "Location contour confidence level" },
	{ PilNone,   "",   "" }
	};


class MultiParams: public PilParams
{
public:
	MultiParams(): PilParams(c_params) {}
	void Print() { PilParams::Print(c_params); }
};


class AppScreen
{
public:
	AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "#### AG_Multi5-" << ROIMULTI_VERSION << " - " << ROIMULTI_DATE <<  " - A.C., T.C., A.T. #####"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	}

	~AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "##########  Task AG_Multi5........... exiting ###################"<< endl;
	cout << "#################################################################"<< endl;
	}
};



int main(int argc, char *argv[])
{
AppScreen appScreen;

MultiParams mPars;
if (!mPars.Load(argc, argv))
	return -1;
cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;
mPars.Print();

MapList maplist;
int mapCount = maplist.Read(mPars["maplist"]);
if (!mapCount) {
	cerr << "File " << mPars.GetStrValue("maplist") << " missing or empty" << endl;
	return -1;
	}

MapData mapData;
if (!mapData.Load(maplist))
	return -1;

SourceDataArray srcArr = ReadSourceFile(mPars["srclist"]);
if (!srcArr.Count())
	cout << "Warning: no point sources loaded" << endl;

// ExpCorr expCorr("None"); /// zzz To remove

RoiMulti roiMulti;
if (!roiMulti.SetPsf(mPars["psdfile"], mPars["sarfile"], mPars["edpfile"]))
	return -1;
if (!roiMulti.SetMaps(mapData , mPars["galmode"], mPars["isomode"]))
	return -1;

const char* outfilename = mPars["outfile"];
string fileName;

fileName = outfilename;
fileName += ".log";

roiMulti.SetLogfile(fileName.c_str());
if (roiMulti.DoFit(srcArr, mPars["ranal"], mPars["ulcl"], mPars["loccl"], 1))
	return -1;

roiMulti.Write(outfilename);
roiMulti.WriteSources(outfilename);
roiMulti.WriteHtml(outfilename);

fileName = outfilename;
fileName += "2";roiMulti.Write(fileName.c_str(), false);

return 0;
}
