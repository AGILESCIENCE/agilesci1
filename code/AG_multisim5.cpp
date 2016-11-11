/*
 * Copyright (c) 2010-2016
 *     Andrew Chen, Tommaso Contessi (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#include <TRandom3.h>
#include <RoiMulti5.h>
#include <PilParams.h>
#include <Eval.h>
#include <FitsUtils.h>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;

const char* startString = {
"#################################################################\n"
"### Task AG_multisim5 v1.2.0 - A.C., T.C., A.B., A.Z.         ###"
};

const char* endString = {
"### Task AG_multisim5 exiting ................................ ###\n"
"##################################################################"
};

const PilDescription paramsDescr[] = {
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
	{ PilString, "resmatrices", "Response matrices" },
	{ PilString, "respath", "Response matrices path" },
	{ PilNone,   "", "" }
};

enum { Concise=1, SkipAnalysis=2, DoubleAnalysis=4, SaveMaps=8 };

AgileMap SumMaps(const AgileMap* mapArr, int offset, int block) {
	AgileMap m(mapArr[offset]);
	for (int i=offset+1; i<offset+block; ++i)
		m += mapArr[i];
	return m;
}

AgileMap SumExposure(const MapMaps& maps, int offset, int block) {
	AgileMap m(maps.ExpMap(offset));
	for (int i=offset+1; i<offset+block; ++i)
		m += maps.ExpMap(i);
	return m;
}

int main(int argc,char **argv) {
	cout << startString << endl;

	PilParams params(paramsDescr);
	if (!params.Load(argc, argv))
		return EXIT_FAILURE;

	cout << endl << "INPUT PARAMETERS:" << endl;
	params.Print();

	int opmode = params["opmode"];
	int block = params["block"];
	int nruns = params["nruns"];
	int seed = params["seed"];
	const char* sarfilename = params["sarfile"];
	const char* edpfilename = params["edpfile"];
	const char* psdfilename = params["psdfile"];
	const char* maplistsimname = params["maplistsim"];
	const char* srclistsim = params["srclistsim"];
	const char* outfilename = params["outfile"];
	const char* maplistanalysisname = params["maplistanalysis"];
	const char* srclistanalysis = params["srclistanalysis"];
	const char* resmatrices = params["resmatrices"];
	const char* respath = params["respath"];
	double ranal = params["ranal"];
	int galmode = params["galmode"];
	int isomode = params["isomode"];
	double ulcl = params["ulcl"];
	double loccl = params["loccl"];

	if (seed)
		SetSeed(seed);

	MapList maplistsim;
	if (!maplistsim.Read(maplistsimname)) {
		cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
		cout << endString << endl;
		return -1;
	}

	MapList maplistanalysis;
	if ((opmode&SkipAnalysis)==0) {
		if (!maplistanalysis.Read(maplistanalysisname)) {
			cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
			cout << endString << endl;
			return -1;
		}
	}

	if (block>0 && (opmode&SkipAnalysis)==0) {
		if (maplistsim.Count()!=maplistanalysis.Count()) {
			cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
			cout << "ERROR: The two map lists must have the same number of rows" << endl;
			cout << endString << endl;
			return -1;
		}
	}
	if (block>maplistsim.Count())
		block = maplistsim.Count();

	MapData mapData;
	if (!mapData.Load(maplistsim, true)) {
		cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
		cout << "Error reading the simulation map list" << endl;
		cout << endString << endl;
		return -1;
	}
	MapData mapDataAna;
	if ((opmode&SkipAnalysis)==0) {
		if (!mapDataAna.Load(maplistanalysis, true)) {
			cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
			cout << "Error reading the analysis map list" << endl;
			cout << endString << endl;
			return -1;
		}
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
			cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
			cout << "ERROR setting PSF data" << endl;
			cout << endString << endl;
			return -1;
		}
#else
	cout << "Multiple RoiMulti instance" << endl;
#endif

	double sumTS = 0;
	int analysisCount = 0;

	for (int i=0; i<nruns; ++i) {
#ifndef _SINGLE_INSTANCE_
		RoiMulti roiMulti;
		if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
			cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
			cout << "ERROR setting PSF data" << endl;
			cout << endString << endl;
			return -1;
		}
#endif
		cout << endl << "AG_Multisim loop #" << i+1 << endl << endl;

		mapData.MapCoeff::Load(maplistsim);
		roiMulti.SetMaps(mapData);
		cout << "New count maps simulation array size=" << mapData.Count() << endl;
		AgileMap* simArr = roiMulti.NewSimulationArray(srcSimArr); // simArr size == maplistsim size

		if (block) {
			int last = mapData.Length()-block;
			cout << "Using block size=" << block << endl;
			for (int j=0; j<=last; ++j) {
				cout << endl << "Summing maps from " << j+1 << " to " << j+block << " [loop " << i+1 << "]" << endl << endl;
				AgileMap ctsMap = SumMaps(simArr, j, block);
				AgileMap expMap = SumExposure(mapData, j, block);
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
					AgileMap gasMap;
					std::stringstream ss;
					ss << respath << "/" << expMap.GetEmin() << "_" << expMap.GetEmax() << "." << resmatrices << ".disp.conv.sky.gz";
                    std::string diffuseFile = ss.str();
					int status = eval::EvalGasMap(gasMap, expMap, diffuseFile.c_str(), diffuseFile.c_str());
					if(status) {
						cout << "Error during gas map evaluation" << endl;
						cout << "AG_multisim5..................... exiting AG_multisim5 ERROR:" << endl;
						fits_report_error(stdout, status);
						return status;
					}
					cout << endl << "AG_Multisim: Analysis step" << endl << endl;
					MapData analysisMaps(ctsMap, expMap, gasMap, 0, 1, 1);
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

	cout << endString << endl;

	return 0;
}

