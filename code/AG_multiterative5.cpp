////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       Alike
//       Release: V 1.0 -  31/Jul/2009
//	 Release: V 1.0 Andrea Bulgarelli, developed by Tomaso Contessi
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       31/Jul/2009
//                      First release: V1.0
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "RoiMulti5.h"
#include "PilParams.h"

using namespace std;


const PilDescription c_params[] = {
	{ PilString, "maplist", "Map list file name" },
	{ PilString, "sarfile", "SAR file name" },
	{ PilString, "edpfile", "EDP file name" },
	{ PilString, "psdfile", "PSD file name" },
//	{ PilString, "expcorrfile", "Exposure Correction file name" },
	{ PilString, "srcscanlist", "Source scan list" },
	{ PilInt,    "scaniterations", "Max Scan iterations" },
	{ PilReal,   "sqrtsthreshold", "Sqrt(TS) threshold for the entire scan list" },
	{ PilReal,   "distthreshold", "Min distance threshold" },
	{ PilReal,   "fixdistance", "Min distance for fixed sources" },
	{ PilReal,   "minsourcesqrts", "Required sqrt(TS) for each scan source fitting" },
	{ PilReal,   "ranal", "Radius of analysis" },
	{ PilInt,    "galmode", "Diffuse emission mode" },
	{ PilInt,    "isomode", "Isotropic emission mode" },
	{ PilString, "srclist", "Sources list" },
	{ PilString, "outfile", "Output file name" },
	{ PilReal,   "ulcl",    "Upper limit confidence level" },
	{ PilReal,   "loccl",   "Location contour confidence level" },
	{ PilInt,    "fixflagscan", "Fixflag for new sources" },
	{ PilInt,    "fixflagstep2", "Fixflag for second step" },
	{ PilNone,   "",   "" }
	};


class MultiIterParams: public PilParams
{
public:
	MultiIterParams(): PilParams(c_params) {}
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
	char* srcfilename,
	char* outfilename,
	char* expcorrfilename,

	char*   srcscanfilename,
	double& sqrTsThreshold,
	int&    maxIterations,

	double& distanceThreshold,
	double& fixdistance,
	double& minSourceSqrTS,

	double& ranal,
	double& ulcl,
	double& loccl,
	int& galmode,
	int& isomode,
	int& fixflagscan,
	int& fixflagstep2)
{
int status = PILInit(argc, argv);
int numpar = 0;
status = status || PILGetNumParameters(&numpar);
GetPilParam("maplist", maplistname, &status);
GetPilParam("sarfile", sarfilename, &status);
GetPilParam("edpfile", edpfilename, &status);
GetPilParam("psdfile", psdfilename, &status);
GetPilParam("expcorrfile", expcorrfilename, &status);

GetPilParam("srcscanlist", srcscanfilename, &status);
GetPilParam("scaniterations", &maxIterations, &status);
GetPilParam("sqrtsthreshold", &sqrTsThreshold, &status);

GetPilParam("distthreshold", &distanceThreshold, &status);
GetPilParam("fixdistance", &fixdistance, &status);
GetPilParam("minsourcesqrts", &minSourceSqrTS, &status);

GetPilParam("ranal", &ranal, &status);
GetPilParam("galmode", &galmode, &status);
GetPilParam("isomode", &isomode, &status);
GetPilParam("srclist", srcfilename, &status);
GetPilParam("outfile", outfilename, &status);
GetPilParam("ulcl", &ulcl, &status);
GetPilParam("loccl", &loccl, &status);
GetPilParam("fixflagscan", &fixflagscan, &status);
GetPilParam("fixflagstep2", &fixflagstep2, &status);

PILClose(status);
return status;
}


static void PrintInput(
	const char* maplistname,
	const char* sarfilename,
	const char* edpfilename,
	const char* psdfilename,
	const char* srcfilename,
	const char* outfilename,
	char* expcorrfilename,

	char*  srcscanfilename,
	double sqrTsThreshold,
	int    maxIterations,
	double distanceThreshold,
	double fixdistance,
	double minSourceSqrTS,

	double ranal,
	double ulcl,
	double loccl,
	int galmode,
	int isomode,
	int fixflagscan,
	int fixflagstep2)
{
cout << endl << "INPUT PARAMETERS:" << endl << endl;
cout << "Map file name = " << maplistname << endl;
cout << "SAR file name : "<< sarfilename << endl;
cout << "EDP file name : "<< edpfilename << endl;
cout << "PSD file name : "<< psdfilename << endl;
cout << "Exposure Correction file name : "<< expcorrfilename << endl;
cout << "Radius of analysis : " << ranal << endl;

cout << "Galactic parameter mode : " <<  galmode;
if (galmode==0)
	cout << " (constant)" << endl;
else if (galmode==1)
	cout << " (from file)" << endl;
else if (galmode==2)
	cout << " (variable)" << endl;
else if (galmode==3)
	cout << " (single variable parameter)" << endl;

cout << "Isotropic parameter mode : " <<  isomode;
if (isomode==0)
	cout << " (constant)" << endl;
else if (isomode==1)
	cout << " (from file)" << endl;
else if (isomode==2)
	cout << " (variable)" << endl;
else if (isomode==3)
	cout << " (single variable parameter)" << endl;

cout << "Upper limit confidence level: " << ulcl << endl;
cout << "Location contour confidence level: " << loccl << endl;
cout << "source list : " <<  srcfilename << endl;	
cout << "Output file name prefix : " <<  outfilename << endl;

cout << "Source Scanlist file name : " <<  srcscanfilename << endl;
cout << "Min sqrt(TS) for the scan list at each iteration : " <<  sqrTsThreshold << endl;
cout << "Max Scan iterations : " <<  maxIterations << endl;
cout << "Distance max treshold for scanning : " <<  distanceThreshold << endl;
cout << "Max distance for fixing the source position : " <<  fixdistance << " degrees" << endl;
cout << "Min sqrt(TS) for each scan source fitting : " <<  minSourceSqrTS << endl;
cout << "Fixflag for new sources : " <<  fixflagscan << endl;
cout << "Fixflag for second step : " <<  fixflagstep2 << endl;
cout << endl << endl;
}
*/

static string CycleNumber(int n)
{
char buffer[16]; 
sprintf(buffer, "_%02d", n);
return string(buffer);
}


#define AlikeSphdistDeg SphDistDeg

static double MinDistance(const SourceDataArray& baseSrcArr, const SourceData& tryData)
{
double minDistance = 100000;
int count = baseSrcArr.Count();
for (int i=0; i<count; ++i) {
	double distance = AlikeSphdistDeg(baseSrcArr[i].srcL, baseSrcArr[i].srcB, tryData.srcL, tryData.srcB);
	if (distance<minDistance)
		minDistance = distance;
	}
return minDistance;
}


static void ResetDistantFlags(SourceDataArray& srcArr, double srcL, double srcB, double maxDistance)
{
int count = srcArr.Count();
for (int i=0; i<count; ++i)
	if (AlikeSphdistDeg(srcArr[i].srcL, srcArr[i].srcB, srcL, srcB)>maxDistance)
		srcArr[i].fixflag = 0;
}

static void RestoreOriginalFlags(SourceDataArray& changedArr, const SourceDataArray& originalArr)
{
int count = originalArr.Count();
for (int i=0; i<count; ++i)
	changedArr[originalArr[i].label].fixflag = originalArr[i].fixflag;
}


static void PrintSourceSummary(ofstream& oFile, const SourceDataArray& srcArr)
{
int srcCount = srcArr.Count();
for (int i=0; i<srcCount; ++i) {
	const SourceData& srcData = srcArr[i];
	oFile << srcData.label << "\t";
	oFile << srcData.srcL << "\t";
	oFile << srcData.srcB << "\t";
	oFile << srcData.TS << "\t";
	oFile << srcData.flux << "\t";
	oFile << srcData.index << "\t";
	oFile << srcData.fixflag << "\t";
	oFile << srcData.minTS << "\t";
	oFile << srcData.gal << "\t";
	oFile << srcData.iso << "\t";
	oFile <<  endl;
	}
}


static void ScrPrint(const char* msg, int num)
{
cout << "AG_Multiterative4: " << msg << " " << num << endl;
}


class AppScreen
{
public:
	AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "#### AG_multiterative4 v.1.1 - 15/03/2011 - A.B., T.C., A.C. ####"<< endl;
	cout << "#################################################################"<< endl;
	cout << "#################################################################"<< endl;
	}

	~AppScreen()
	{
	cout << "#################################################################"<< endl;
	cout << "#######  Task AG_Multiterative4........ exiting #################"<< endl;
	cout << "#################################################################"<< endl;
	}
};


#define AlikeMap AgileMap



int main(int argc, char *argv[])
{
AppScreen appScreen;

MultiIterParams mPars;
if (!mPars.Load(argc, argv))
	return -1;
cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;
mPars.Print();

const char* maplistname = mPars["maplist"];
const char* sarfilename = mPars["sarfile"];
const char* edpfilename = mPars["edpfile"];
const char* psdfilename = mPars["psdfile"];
//const char* expcorrfilename = mPars["expcorrfile"];

const char* srcscanfilename = mPars["srcscanlist"];
int maxIterations = mPars["scaniterations"];

double sqrTsThreshold = mPars["sqrtsthreshold"];
double distanceThreshold = mPars["distthreshold"];
double fixdistance = mPars["fixdistance"];
double minSourceSqrTS = mPars["minsourcesqrts"];

double ranal = mPars["ranal"];
int galmode = mPars["galmode"];
int isomode = mPars["isomode"];

const char* srcfilename = mPars["srclist"];
const char* outfilename = mPars["outfile"];

double ulcl = mPars["ulcl"];
double loccl = mPars["loccl"];

FixFlag fixflagscan = mPars["fixflagscan"];
FixFlag fixflagstep2 = mPars["fixflagstep2"];





MapList maplist;
int mapCount = maplist.Read(maplistname);
if (!mapCount) {
	cerr << "ERROR: File " << maplistname << " missing or empty" << endl;
	return -1;
	}
MapData mapData;
if (!mapData.Load(maplist))
	return -1;


double tsThreshold = sqrTsThreshold*sqrTsThreshold;
double minSourceTS = minSourceSqrTS*minSourceSqrTS;

// ExpCorr expCorr(expcorrfilename);

SourceDataArray baseSrcArr = ReadSourceFile(srcfilename);
if (!baseSrcArr.Count())
	cerr << "Warning: File " << srcfilename << " not found or empty" << endl;
SourceDataArray scanSrcArr = ReadSourceFile(srcscanfilename);
if (!scanSrcArr.Count()) {
	cerr << "ERROR: File " << srcscanfilename << " not found or empty" << endl;
	return -1;
	}
int tryCount = scanSrcArr.Count();
int* originalFlags = new int[baseSrcArr.Count()+maxIterations];
double galc = 0;
double isoc = 0;

RoiMulti roiMulti;
if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
	cerr << "ERROR accessing PSF related files" << endl;
	return -1;
	}
if (!roiMulti.SetMaps(mapData, galmode, isomode)) {
	cerr << "ERROR accessing the map list" << endl;
	return -1;
	}

AlikeMap modelMap(maplist.CtsName(0)); /// Model for flux, ts, and index maps

for (int cycle=0; cycle<maxIterations; ++cycle) {
	//Log file
	string outlogname(outfilename);
	outlogname += ".log";
	outlogname += CycleNumber(cycle);
	ofstream logFile(outlogname.c_str());

	ScrPrint("Multiterative cycle #", cycle);
	logFile << "!! Start Step One of cycle " << cycle << endl;

	/// STEP ONE: try adding each source of scanSrcArr to baseSrcArr and find the best
	SourceData bestTry;

	int baseCount = baseSrcArr.Count();
	for (int i=0; i<baseCount; ++i) {
		originalFlags[i] = baseSrcArr[i].fixflag;
		//baseSrcArr[i].fixflag = 1;		/// All the base sources will have fixflag=1 in the step one
		}

	//fixflag = 1 is used in step one for the source of the scan list under test
	//For step one the other source has
	// - the original fixflag for the input list sources
	// - the original fixflag for the scan list sources
	// - but if the source are too far with respect to the try source, the fixflag = 0. This means that it is necessary that the flux is present for the sources of the input list
	//fixflag = 2 is used in step two for the best source of the scan list
	//For step two the other sources the rules are the same of step one
	double maxTS = 0;
	int tryIndexNum = 1;
// double tryIndex[3] = {2.2, 1.8, 1.5};
	double tryIndex[1] = {2.1};
		
	for (int i=0; i<tryCount; ++i) {
		SourceData tryData = scanSrcArr[i];

		logFile << "================================" << endl;
		logFile << "Analyze source " << i << " (" << tryData.label << ", " << tryData.srcL << ", " << tryData.srcB << ", " << tryData.TS << ", " << tryData.flux << ") " << endl;

		SourceData bestCurrentTry;
		double maxTSCurrentTry = -999;
		double minTSScan = 0.0;
		double minDistance = MinDistance(baseSrcArr, tryData);

		// Make a try if the source is not too close and if its previous TS was big enough
		// Be carefull, tryData.TS > X is applied also in the case of scan spectral index.
		// This means that if the source with starting spectral index has TS < X, no other
		// spectral index are applied

		if (minDistance>distanceThreshold) {
			//La sorgente è abbastanza lontano dalle altre, valuta l'index scan
			for (int j=0; j<tryIndexNum; j++) {				
				tryData.index = tryIndex[j];
				
				//almeno nel primo ciclo vanno valutate tutte, perche' tutte con TS=0
				if(cycle==0 || tryData.TS >= minTSScan) {
					/// Make a copy of the base sources setting fixflag=0 for those too far apart
					logFile << "* Starting with index ... " << tryData.index << " and minDistance " << minDistance << " and current maxTS " << maxTS << endl;
					
					ScrPrint("Evalueate source with spectral index #", tryData.index);
					
					SourceDataArray tryArr(baseSrcArr);
					ResetDistantFlags(tryArr, tryData.srcL, tryData.srcB, fixdistance);
	
					/// Assign a new name and fixflag for this try
					string tryName = tryData.label + CycleNumber(cycle);
					tryData.label = tryName;
					tryData.fixflag = fixflagscan; //test the current position
					tryData.flux = 0; //questo perche' se non fa lo step 2 non calcola nemmeno l'UL
					
					tryArr.Append(tryData);
					tryArr.Print(logFile);
					/// tryArr.Print();	/// zzz debug
					roiMulti.DoFit(tryArr, ranal, ulcl, loccl, 0, tryData.label.c_str(), minSourceTS);
					galc = roiMulti.GetGalactic(0).GetCoeff();
					isoc = roiMulti.GetIsotropic(0).GetCoeff();
					tryArr = roiMulti.GetFitData();
					/// tryArr.Print(true);	/// zzz debug
					tryData = tryArr[tryName];	/// Get the data for the current try
					if(tryData.TS < 0) {
						tryData.TS = 0;
					}
					logFile << "Result: " << i << " (" << tryData.label << ", " << tryData.srcL << ", " << tryData.srcB << ", " << tryData.TS << ", " << tryData.flux << ") with index " << tryData.index << endl;

					tryData.fixflag = scanSrcArr[i].fixflag;	/// Restore the original fixlag
					tryData.gal= galc;
					tryData.iso= isoc;
					//select a local best try
					if (tryData.TS > maxTSCurrentTry) {
						maxTSCurrentTry = tryData.TS;
						bestCurrentTry = tryData;
					}
					//select the global best try
					if (tryData.TS>maxTS) {
						maxTS = tryData.TS;
						bestTry = tryData;
						logFile << "################ Found a new TS max " << maxTS << endl;
					}
					tryData.label = scanSrcArr[i].label;	/// Restore the original name
					bestCurrentTry.label = scanSrcArr[i].label; //rimetto a posto la label originale
					scanSrcArr[i] = bestCurrentTry;
					//non effettuare la ricerca per altri indici, cambia sorgente
					if(tryData.TS < 1)
						break;
				} else {
					//in ogni caso, anche se non si rivaluta per un TS troppo piccolo,
					//il suo TS va registrato e valutato comunque se e' un max
					if (tryData.TS>maxTS) {
						maxTS = tryData.TS;
						bestTry = tryData;
						logFile << "################" << endl;
						logFile << "Found a new TS max without direct evaluation of the source " << maxTS << endl;
					}
				}	
			}
		} else {
			//la sorgente e' troppo vicina ad un'altra. Non valutarla,
			//E poi che faccio? Tengo il vecchio TS? Perche' potrebbe essere
			//rivalutabile nei cicli successivi.	
			logFile << "Source too near to another source: " << minDistance << endl;
		}
		logFile << "Result for this index scan: " << i << " (" << scanSrcArr[i].label << ", " << scanSrcArr[i].srcL << ", " << scanSrcArr[i].srcB << ", " << scanSrcArr[i].TS << ", " << scanSrcArr[i].flux << ") with index " << scanSrcArr[i].index << endl;
	}
	for (int i=0; i<baseCount; ++i)
		baseSrcArr[i].fixflag = originalFlags[i];	/// Restore all the original flags to the base sources

	///AB, salva lo scan src array dello step corrente
	string outfname(outfilename);
	outfname += ".lst";
	outfname += CycleNumber(cycle);

	ofstream asciiFile(outfname.c_str(), ios::app);
	if(asciiFile.is_open())
		PrintSourceSummary(asciiFile, scanSrcArr);
	asciiFile.close();

	logFile << "!! Start saving output files of cycle " << cycle << endl;

	/// Print TS and flux maps
	{
		AlikeMap fluxMap(modelMap);
		AlikeMap tsMap(modelMap);
		AlikeMap siMap(modelMap);
		AlikeMap gasMap(modelMap);
		AlikeMap isoMap(modelMap);
		fluxMap.Zero();
		tsMap.Zero();
		siMap.Zero();
		gasMap.Zero();
		isoMap.Zero();
		for (int i=0; i<tryCount; ++i) {
			SourceData tryData = scanSrcArr[i];
			int row, col;
			bool inside = fluxMap.GetRowCol(tryData.srcL, tryData.srcB, &row, &col);
			if (inside) {
				fluxMap(row, col) = (Double_t)tryData.flux; //AB1
				tsMap(row, col) = (Double_t)tryData.TS; //AB1
				siMap(row, col) = (Double_t)tryData.index; //AB1
				gasMap(row, col) = tryData.gal;
				isoMap(row, col) = tryData.iso;
				}
			}
		string outfname2(outfname);
		outfname2 += ".flux.fits.gz";
		fluxMap.Write(outfname2.c_str());
		string outfname3(outfname);
		outfname3 += ".TS.fits.gz";
		tsMap.Write(outfname3.c_str()); //AB1
		string outfname4(outfname);
		outfname4 += ".SI.fits.gz";
		siMap.Write(outfname4.c_str()); //AB1
	
		string outfname5(outfname);
		outfname5 += ".GAS.fits.gz";
		gasMap.Write(outfname5.c_str()); //AB1
		string outfname6(outfname);
		outfname6 += ".ISO.fits.gz";
		isoMap.Write(outfname6.c_str()); //AB1
		}

	logFile << "!! Start Step Two of cycle " << cycle << endl;
	/// STEP TWO: Add the best try found to the base sources with all the original flags
	if (maxTS>=tsThreshold) {
		SourceDataArray tryArr(baseSrcArr);
		SourceData bestTryOriginal = bestTry;
		ResetDistantFlags(tryArr, bestTry.srcL, bestTry.srcB, fixdistance);

		string tryName = bestTry.label;
		int originalFixflag = bestTry.fixflag;
		double startL = bestTry.srcL;
		double startB = bestTry.srcB;
		bestTry.fixflag = fixflagstep2;
		tryArr.Append(bestTry);

		/// Printing the input
		cout << "Input sources for Step Two:" << endl;
		tryArr.Print(cout);
		logFile << "Input sources for Step Two:" << endl;
		tryArr.Print(logFile);
		
		string inputStepFname(outfilename);
		inputStepFname += ".inputstep2";
		inputStepFname += CycleNumber(cycle);
		ofstream inputStepFile(inputStepFname.c_str());
		tryArr.Print(inputStepFile);

		string fileName(outfilename);
		fileName += CycleNumber(cycle);
		roiMulti.DoFit(tryArr, ranal, ulcl, loccl, 1);
		tryArr = roiMulti.GetFitData();
		galc = roiMulti.GetGalactic(0).GetCoeff();
		isoc = roiMulti.GetIsotropic(0).GetCoeff();
		roiMulti.Write(fileName.c_str());
		roiMulti.WriteSources(fileName.c_str(), true);
		roiMulti.WriteHtml(fileName.c_str());

		/// Check if the position calculated of the new source is too far.
		double distanceOldNewPosition = AlikeSphdistDeg(tryArr[tryName].srcL, tryArr[tryName].srcB, startL, startB);
		if(distanceOldNewPosition > 4.0) {
			cout << "WARNING: New position found is too far, try with fixflag = 1" << endl;
			logFile << "WARNING: New position found is too far, try with fixflag = 1" << endl;
			tryArr[tryName]  = bestTryOriginal;
			tryArr[tryName].fixflag = 1;
			tryArr.Print(cout);
			tryArr.Print(logFile);
			roiMulti.DoFit(tryArr, ranal, ulcl, loccl, 1);
			tryArr = roiMulti.GetFitData();
			galc = roiMulti.GetGalactic(0).GetCoeff();
			isoc = roiMulti.GetIsotropic(0).GetCoeff();
			roiMulti.Write(fileName.c_str());
			roiMulti.WriteSources(fileName.c_str(), true);
			roiMulti.WriteHtml(fileName.c_str());
			}
		/// Restore the fixflags
		tryArr[tryName].fixflag = originalFixflag;
		tryArr[tryName].gal= galc;
		tryArr[tryName].iso= isoc;
		RestoreOriginalFlags(tryArr, baseSrcArr);
		baseSrcArr = tryArr;
		}
	else
		break;
	}
delete[] originalFlags;
return 0;
}
