////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_multisim
//       Release: BUILD 21 -  30/Dec/2010
//       Contributors: 
//       Author: Andrew Chen, Tomaso Contessi (IASF-Milano)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       30/Dec/2010
//                      First release: V0.1
//       		Author: Andrew Chen, Tomaso Contessi (IASF-Milano)
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////



#include <TRandom3.h>

#include "RoiMulti5.h"
#include "PilParams.h"


using namespace std;


const PilDescription c_params[] = {
	{ PilInt,    "opmode", "Operation Mode" },
	{ PilInt,    "block", "Block" },
	{ PilInt,    "nruns", "Number of runs" },
	{ PilInt,    "seed", "Seed" },
	{ PilString, "sarfile", "SAR file name" },
	{ PilString, "edpfile", "EDP file name" },
	{ PilString, "psdfile", "PSD file name" },


	{ PilString, "maplistsim", "Map list for simulation" },
	{ PilString, "srclistsim", "Source list for simulation" },
	{ PilString, "outfile", "Output file name" },

	{ PilString, "maplistanalysis", "Map list for analysis" },
	{ PilString, "srclistanalysis", "Source list for analysis" },
	{ PilReal,   "ranal", "Radius of analysis" },
	{ PilInt,    "galmode", "Diffuse emission mode" },
	{ PilInt,    "isomode", "Isotropic emission mode" },
	{ PilReal,   "ulcl",    "Upper limit confidence level" },
	{ PilReal,   "loccl",   "Location contour confidence level" },
	{ PilNone,   "",   "" }
	};


class MultiSimParams: public PilParams
{
public:
	MultiSimParams(): PilParams(c_params) {}
	void Print() { PilParams::Print(c_params); }
};


enum { Concise=1, SkipAnalysis=2, DoubleAnalysis=4, SaveMaps=8 };


/**
static int DoPIL(
	int argc,
	char *argv[],

	int& opmode,
	int& block,
	int& nruns,
	int& seed,

	char* sarfilename,
	char* edpfilename,
	char* psdfilename,

	char* maplistsimname,
	char* srclistsim,
	char* outfilename,

	char* maplistanalysisname,
	char* srclistanalysis,
	double& ranal,
	int& galmode,
	int& isomode,
	double& ulcl,
	double& loccl)
{
int status = PILInit(argc, argv);
int numpar = 0;
status = status || PILGetNumParameters(&numpar);

GetPilParam("opmode", &opmode, &status);
GetPilParam("block", &block, &status);
if (block<0)
	block = -block;
GetPilParam("nruns", &nruns, &status);
GetPilParam("seed", &seed, &status);

GetPilParam("sarfile", sarfilename, &status);
GetPilParam("edpfile", edpfilename, &status);
GetPilParam("psdfile", psdfilename, &status);

GetPilParam("maplistsim", maplistsimname, &status);
GetPilParam("srclistsim", srclistsim, &status);

GetPilParam("outfile", outfilename, &status);

if ((opmode&SkipAnalysis)==0) {
	GetPilParam("maplistanalysis", maplistanalysisname, &status);
	GetPilParam("srclistanalysis", srclistanalysis, &status);
	GetPilParam("ranal", &ranal, &status);
	GetPilParam("galmode", &galmode, &status);
	GetPilParam("isomode", &isomode, &status);
	GetPilParam("ulcl", &ulcl, &status);
	GetPilParam("loccl", &loccl, &status);
	}

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
	int opmode,
	int block,
	int nruns,
	int seed,

	const char* sarfilename,
	const char* edpfilename,
	const char* psdfilename,

	const char* maplistsimname,
	const char* srclistsim,

	const char* outfilename,

	const char* maplistanalysisname,
	const char* srclistanalysis,
	double ranal,
	int galmode,
	int isomode,
	double ulcl,
	double loccl)
{
cout << endl << endl << "INPUT PARAMETERS:"<< endl << endl;

cout << "Operation mode: " << opmode;
if (opmode & Concise)
	cout << " Concise output";
if (opmode & SkipAnalysis)
	cout << " No analysis performed";
else if (opmode & DoubleAnalysis)
	cout << " Double Analysis performed";
if (opmode & SaveMaps)
	cout << " Intermediate maps saved";
if (opmode==0)
	cout << " Verbose output";
cout << endl;

cout << "Block: " << block;
if (block==0)
	cout << " (No Block iteration)" << endl;
else
	cout << " (Block iteration)" << endl;

cout << "Number of runs: " << nruns << endl;

if (seed)
	cout << "Seed: " << seed << endl;
else
	cout << "Seed: Random (" << GetSeed() << ")" << endl;

cout << endl << "Calibration" << endl;
cout << "SAR file name: " << sarfilename << endl;
cout << "EDP file name: " << edpfilename << endl;
cout << "PSD file name: " << psdfilename << endl;

cout << endl << "Simulation" << endl;
cout << "Simulation Map list: " << maplistsimname << endl;
cout << "Simulation Source list: " <<  srclistsim << endl;
cout << "Output file prefix: " <<  outfilename << endl;

cout << endl << "Analysis";
if (opmode&SkipAnalysis)
	cout << " not performed" << endl;
else {
	cout << endl;
	cout << "Analysis Map list: " << maplistanalysisname << endl;
	cout << "Analysis Source list: " <<  srclistanalysis << endl;	
	cout << "Radius of analysis: " << ranal << endl;
	cout << "Galactic parameter mode: "; PrintDiffMode(galmode);
	cout << "Isotropic parameter mode: "; PrintDiffMode(isomode);
	cout << "Upper limit confidence level: " << ulcl << endl;
	cout << "Location contour confidence level: " << loccl << endl;
	}
cout << endl;
}
*/


class AppScreen
{
public:
	AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "###### AG_multisim4  v.1.0 - 25/11/2011 - A.C., T.C. ############"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	cout << endl << "RoiMulti version: " << ROIMULTI_VERSION << " - " << ROIMULTI_DATE<< endl;
	}

	~AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "##########  Task AG_multisim4.......... exiting #################"<< endl;
	cout << "#################################################################"<< endl;
	}
};


#define AlikeMap AgileMap
/**
AlikeMap SumMaps(const AlikeMap* mapArr, int offset, int count)
{
AlikeMap m(mapArr[offset]);
for (int i=offset+1; i<offset+count; ++i)
	m.TMatrixD::operator+=(mapArr[i]);
return m;
}

AlikeMap SumExposure(const MapMaps& maps, int offset, int count)
{
AlikeMap m(maps.ExpMap(offset));
for (int i=offset+1; i<offset+count; ++i)
	m.TMatrixD::operator+=(maps.ExpMap(i));
return m;
}
*/
AgileMap SumMaps(const AgileMap* mapArr, int offset, int count)
{
AgileMap m(mapArr[offset]);
for (int i=offset+1; i<offset+count; ++i)
	m += mapArr[i];
return m;
}

AgileMap SumExposure(const MapMaps& maps, int offset, int count)
{
AgileMap m(maps.ExpMap(offset));
for (int i=offset+1; i<offset+count; ++i)
	m += maps.ExpMap(i);
return m;
}


int main(int argc,char **argv)
{
AppScreen appScreen;


MultiSimParams mPars;
if (!mPars.Load(argc, argv))
	return -1;
cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;
mPars.Print();

int opmode = mPars["opmode"];
int block = mPars["block"];
int nruns = mPars["nruns"];
int seed = mPars["seed"];
const char* sarfilename = mPars["sarfile"];
const char* edpfilename = mPars["edpfile"];
const char* psdfilename = mPars["psdfile"];
const char* maplistsimname = mPars["maplistsim"];
const char* srclistsim = mPars["srclistsim"];
const char* outfilename = mPars["outfile"];
const char* maplistanalysisname = mPars["maplistanalysis"];
const char* srclistanalysis = mPars["srclistanalysis"];



double ranal = mPars["ranal"];
int galmode = mPars["galmode"];
int isomode = mPars["isomode"];
double ulcl = mPars["ulcl"];
double loccl = mPars["loccl"];


	if (seed)
		SetSeed(seed);

	MapList maplistsim;
	if (!maplistsim.Read(maplistsimname))
		return -1;
	MapList maplistanalysis;
	if ((opmode&SkipAnalysis)==0)
		if (!maplistanalysis.Read(maplistanalysisname))
			return -1;

	if (block>0 && (opmode&SkipAnalysis)==0)
		if (maplistsim.Count()!=maplistanalysis.Count()) {
			cerr << "ERROR: The two map lists must have the same number of rows" << endl;
			return -1;
			}
	if (block>maplistsim.Count())
		block = maplistsim.Count();

	MapData mapData;
	if (!mapData.Load(maplistsim, true)) {
		cerr << "Error reading the simulation map list" << endl;
		return -1;
		}
	MapData mapDataAna;
	if ((opmode&SkipAnalysis)==0)
		if (!mapDataAna.Load(maplistanalysis, true)) {
			cerr << "Error reading the analysis map list" << endl;
			return -1;
			}

	SourceDataArray srcSimArr = ReadSourceFile(srclistsim);
	SourceDataArray srcAnaArr;
	if ((opmode&SkipAnalysis)==0)
		srcAnaArr = ReadSourceFile(srclistanalysis);

/// #define _SINGLE_INSTANCE_
#ifdef _SINGLE_INSTANCE_
cout << "Single RoiMulti instance" << endl;
	RoiMulti roiMulti;
	if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
		cerr << "ERROR setting PSF data" << endl;
		return -1;
		}
#else
cout << "Multiple RoiMulti instance" << endl;
#endif

	double sumTS = 0;
	int    analysisCount = 0;

	for (int i=0; i<nruns; ++i) {

#ifndef _SINGLE_INSTANCE_
		RoiMulti roiMulti;
		if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
			cerr << "ERROR setting PSF data" << endl;
			return -1;
			}
#endif
		cout << endl << "AG_Multisim loop #" << i+1 << endl << endl;

		mapData.MapCoeff::Load(maplistsim);
		roiMulti.SetMaps(mapData);
		AlikeMap* simArr = roiMulti.NewSimulationArray(srcSimArr);

		if (block) {
			int last = mapData.Length()-block;
			for (int j=0; j<=last; ++j) {
				cout << endl << "AG_Multisim summing maps from " << j+1 << " to " << j+block << " [loop " << i+1 << "]" << endl << endl;
				AlikeMap ctsMap = SumMaps(simArr, j, block);
				AlikeMap expMap = SumExposure(mapData, j, block);
				/// Writing cts and exp maps
				if (opmode & SaveMaps) {
					char mapName[256];
					sprintf(mapName, "%010d_%03d_%s.cts.gz", i+1, j+1, outfilename);
					if (ctsMap.Write(mapName))
						cerr << "Error writing simulated counts map " << mapName << endl;
					else
						cerr << mapName << " written" << endl;
					sprintf(mapName, "%010d_%03d_%s.exp.gz", i+1, j+1, outfilename);
					if (expMap.Write(mapName))
						cerr << "Error writing simulated counts map " << mapName << endl;
					else
						cerr << mapName << " written" << endl;
					}

				if ((opmode&SkipAnalysis)==0) {
					cout << endl << "AG_Multisim: Analysis step" << endl << endl;
					MapData analysisMaps(ctsMap, expMap, mapDataAna.GasMap(0), 0, 1, 1);
					analysisMaps.MapCoeff::Load(maplistanalysis);
					roiMulti.SetMaps(analysisMaps, galmode, isomode);
					roiMulti.DoFit(srcAnaArr, ranal, ulcl, loccl);
					if (opmode & DoubleAnalysis) {
						cout << endl << "AG_Multisim: Second analysis step" << endl << endl;
						SourceDataArray newSources = roiMulti.GetFitData();
						roiMulti.DoFit(newSources, ranal, ulcl, loccl);
						}

					++analysisCount;
					SourceDataArray results = roiMulti.GetFitData();
					int dataCount = results.Count();
					for (int s=0; s<dataCount; ++s)
						if (results[s].fixflag!=AllFixed)
							sumTS += results[s].TS;

					if (opmode & Concise)
						roiMulti.LogSources(outfilename, i+1);
					else {
						char fileName[256];
						sprintf(fileName, "%010d_%03d_%s", i+1, j+1, outfilename);
						roiMulti.Write(fileName);
						roiMulti.WriteSources(fileName, true, true);
						}
					}
				}
			}
		else {

			if (opmode & SaveMaps) {
				char ctsName[256];
				for (int m=0; m<mapData.Length(); ++m) {
					sprintf(ctsName, "%010d_%03d_%s.cts.gz", i+1, m+1, outfilename);
					if (simArr[0].Write(ctsName))
						cerr << "Error writing simulated counts map " << ctsName << endl;
					else
						cerr << ctsName << " written" << endl;
					}
				}

			if ((opmode&SkipAnalysis)==0) {
				cout << endl << "AG_Multisim: Analysis step" << endl << endl;
				mapDataAna.ReplaceCtsMaps(simArr);

				roiMulti.SetMaps(mapDataAna, galmode, isomode);
				roiMulti.DoFit(srcAnaArr, ranal, ulcl, loccl);

				if (opmode & DoubleAnalysis) {
					cout << endl << "AG_Multisim: Second analysis step" << endl << endl;
					SourceDataArray newSources = roiMulti.GetFitData();
					roiMulti.DoFit(newSources, ranal, ulcl, loccl);
					}

				++analysisCount;
				SourceDataArray results = roiMulti.GetFitData();
				int dataCount = results.Count();
				for (int s=0; s<dataCount; ++s)
					if (results[s].fixflag!=AllFixed)
						sumTS += results[s].TS;

				if (opmode & Concise)
					roiMulti.LogSources(outfilename, i+1);
				else {
					char fileName[256];
					sprintf(fileName, "%010d_%s", i+1, outfilename);
					roiMulti.Write(fileName);
					roiMulti.WriteSources(fileName, true, true);
					}
				}
			}
		delete[] simArr;
		}
	if (analysisCount)
		cout << endl << "AG_Multisim performed the analysis " << analysisCount << " times" << endl << "Total TS: " << sumTS << ", average: " << sumTS/analysisCount << endl << endl;


return 0;
}

