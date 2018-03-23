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
	{ PilInt,    "galmode2", "Diffuse emission optimisation for Loop2" },
	{ PilInt,    "galmode2fit", "Diffuse emission optimisation for Loop2 - fit" },
	{ PilInt,    "isomode2", "Isotropic emission optimisation for Loop2" },
	{ PilInt,    "isomode2fit", "Isotropic emission optimisation for Loop2 - fit" },
	{ PilInt,    "edpcorrection", "EDP cleaning correction" },
	{ PilInt,    "fluxcorrection", "Flux calculation correction for spectral shape" },
	{ PilString, "minimizertype", "Minimizer type" },
	{ PilString, "minimizeralg", "Minimizer algorithm" },
	{ PilInt,    "minimizerdefstrategy", "Minimizer default strategy" },
	{ PilReal,   "mindefaulttolerance", "Minimizer default tolerance"},
	{ PilBool, "expratioevaluation","If 'yes' (or 'y') the exp-ratio evaluation will be enabled."},
	{ PilBool, "isExpMapNormalized","If 'yes' (or 'y') you assert that the exp-map is already normalized. Insert 'no' (or 'n') instead and the map will be normalized before carrying out the exp-ratio evaluation."},
	{ PilReal, "minThreshold", "The lower bound for the threshold level in exp-ratio evaluation"},
	{ PilReal, "maxThreshold", "The upper bound for the threshold level in exp-ratio evaluation"},
	{ PilReal, "squareSize", "The edge degree dimension of the exp-ratio evaluation area"},
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
	cout << "#### AG_Multi6 v2.0.0 - A.C., T.C., A.T., A.B               #####"<< endl;
	cout << "#################################################################"<< endl;
	}

	~AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "##########  Task AG_Multi6........... exiting ###################"<< endl;
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
roiMulti.SetMinimizer(mPars["minimizertype"], mPars["minimizeralg"], mPars["minimizerdefstrategy"], mPars["mindefaulttolerance"]);
roiMulti.SetCorrections(mPars["galmode2"], mPars["galmode2fit"], mPars["isomode2"], mPars["isomode2fit"], mPars["edpcorrection"], mPars["fluxcorrection"]);
if (roiMulti.DoFit(srcArr, mPars["ranal"], mPars["ulcl"], mPars["loccl"], 1))
	return -1;

roiMulti.Write(outfilename);
roiMulti.WriteSources(outfilename, mPars["expratioevaluation"], mPars["isExpMapNormalized"], mPars["minThreshold"], mPars["maxThreshold"], mPars["squareSize"]);
roiMulti.WriteHtml(outfilename, mPars["expratioevaluation"], mPars["isExpMapNormalized"], mPars["minThreshold"], mPars["maxThreshold"], mPars["squareSize"]);

fileName = outfilename;
fileName += ".log";roiMulti.Write(fileName.c_str(), false);

return 0;
}
