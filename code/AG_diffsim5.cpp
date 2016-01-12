////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_diffsim
//       Contributors:
//       Author: Andrea Zoli, Andrea Bulgarelli (IASF-Bologna)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

// TODO3: test with sources

#include <TRandom3.h>
#include <RoiMulti5.h>
#include <PilParams.h>

using namespace std;

const PilDescription c_params[] = {
    { PilInt,    "opmode", "Operation Mode" },
    { PilInt,    "nruns", "Number of runs" },
    { PilInt,    "seed", "Seed" },
    { PilInt,    "radius", "Radius from the map center" },
    { PilString, "sarfile", "SAR file name" },
    { PilString, "edpfile", "EDP file name" },
    { PilString, "psdfile", "PSD file name" },
    { PilString, "maplistsim3", "Map list 3 for simulation" },
    { PilString, "maplistsim7", "Map list 7 for simulation" },
    { PilString, "srclistsim", "Source list for simulation" },
    { PilString, "outfile", "Output file name" },
    { PilNone,   "", "" }
};

class DiffSimParams: public PilParams {
public:
    DiffSimParams() : PilParams(c_params) {}
    void Print() { PilParams::Print(c_params); }
};


enum { Analysis=1, SaveMaps=2 };

class AppScreen {
public:
    AppScreen() {
        cout << "#################################################################"<< endl;
        cout << "######## AG_diffsim  v.1.0 - 11/01/2016 - A.Z., A.B. ############"<< endl;
        cout << "#################################################################"<< endl;
        cout << "#################################################################"<< endl;
        cout << endl << "RoiMulti version: " << ROIMULTI_VERSION << " - " << ROIMULTI_DATE << endl;
    }

    ~AppScreen() {
        cout << "#################################################################"<< endl;
        cout << "##########  Task AG_diffsim............ exiting #################"<< endl;
        cout << "#################################################################"<< endl;
    }
};

int main(int argc,char **argv) {
    AppScreen appScreen;

    DiffSimParams mPars;
    if (!mPars.Load(argc, argv))
        return -1;
    cout << endl << endl << "INPUT PARAMETERS:" << endl << endl;
    mPars.Print();

    int opmode = mPars["opmode"];
    int nruns = mPars["nruns"];
    int seed = mPars["seed"];
    int radius = mPars["radius"];
    const char* sarfilename = mPars["sarfile"];
    const char* edpfilename = mPars["edpfile"];
    const char* psdfilename = mPars["psdfile"];
    const char* maplistsimname3 = mPars["maplistsim3"];
    const char* maplistsimname7 = mPars["maplistsim7"];
    const char* srclistsim = mPars["srclistsim"];
    const char* outfilename = mPars["outfile"];

    if (seed)
        SetSeed(seed);

    bool bothMaps = true;
    if (string(maplistsimname7).compare("None") == 0 ||
        string(maplistsimname7).compare("none") == 0)
        bothMaps = false;

    MapList maplistsim7;
    if (bothMaps && !maplistsim7.Read(maplistsimname7)) {
         return -1;
    }

    MapList maplistsim3;
    if (!maplistsim3.Read(maplistsimname3))
        return -1;

    MapData mapData7;
    if(bothMaps && !mapData7.Load(maplistsim7, true)) {
        cerr << "Error reading the simulation map list 7" << endl;
        return -1;
    }

    MapData mapData3;
    if (!mapData3.Load(maplistsim3, true)) {
        cerr << "Error reading the simulation map list 3" << endl;
        return -1;
    }

    SourceDataArray srcSimArr = ReadSourceFile(srclistsim);

    char fName[256];
    int sizeX = 0;
    int sizeY = 0;
    double centerX = 0.;
    double centerY = 0.;
    ofstream os;
    for (int i=0; i<nruns; ++i) {
        RoiMulti roiMulti;
        if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
            cerr << "ERROR setting PSF data" << endl;
            return -1;
        }

        cout << endl << "AG_diffsim loop #" << i+1 << endl << endl;

        AgileMap* simArr7 = NULL;
        if(bothMaps) {
            mapData7.MapCoeff::Load(maplistsim7);
            roiMulti.SetMaps(mapData7);
            simArr7 = roiMulti.NewSimulationArray(srcSimArr);
        }

        mapData3.MapCoeff::Load(maplistsim3);
        roiMulti.SetMaps(mapData3);
        AgileMap* simArr3 = roiMulti.NewSimulationArray(srcSimArr);

        if (opmode & SaveMaps) {
            sprintf(fName, "%s3_%010d.cts.gz", outfilename, i+1);
            if (simArr3[0].Write(fName))
                cerr << "Error writing simulated counts map 3 " << fName << endl;
            else
                cerr << fName << " written" << endl;
            if(bothMaps) {
                sprintf(fName, "%s7_%010d.cts.gz", outfilename, i+1);
                if (simArr7[0].Write(fName))
                    cerr << "Error writing simulated counts map 7 " << fName << endl;
                else
                    cerr << fName << " written" << endl;
            }
        }
        if (opmode & Analysis) {
            const AgileMap& ctsMap3 = simArr3[0];
            const AgileMap& expMap3 = mapData3.ExpMap(0);

            // first loop analysis init
            if (i == 0) {
                sizeY = expMap3.Dim(0);
                sizeX = expMap3.Dim(1);
                centerX = sizeX / 2.;
                centerY = sizeY / 2.;
                sprintf(fName, "%s_results.txt", outfilename);
                os.open(fName);
            }

            // compute intensity maps
            AgileMap intMap3;
            intMap3.ResizeTo(sizeY, sizeX);
            for (int y=0; y<sizeY; y++) {
                for (int x=0; x<sizeX; x++) {
                    intMap3(y,x) = ctsMap3(y,x) / expMap3(y,x);
                }
            }

            AgileMap intMap7;
            if(bothMaps) {
                const AgileMap& ctsMap7 = simArr7[0];
                const AgileMap& expMap7 = mapData7.ExpMap(0);
                intMap7.ResizeTo(sizeY, sizeX);
                for (int y=0; y<sizeY; y++) {
                    for (int x=0; x<sizeX; x++) {
                        intMap7(y,x) = ctsMap7(y,x) / expMap7(y,x);
                    }
                }
            }

            if (opmode & SaveMaps) {
                sprintf(fName, "%s3_%010d.int.gz", outfilename, i+1);
                if (intMap3.Write(fName))
                    cerr << "Error writing simulated counts map 3 " << fName << endl;
                else
                    cerr << fName << " written" << endl;
                if(bothMaps) {
                    sprintf(fName, "%s7_%010d.int.gz", outfilename, i+1);
                    if (intMap7.Write(fName))
                        cerr << "Error writing simulated counts map 7 " << fName << endl;
                    else
                        cerr << fName << " written" << endl;
                }
            }

            // absolute diff between intensity maps
            AgileMap diffMap;
            diffMap.ResizeTo(sizeY, sizeX);
            if(bothMaps) {
                for (int y=0; y<sizeY; y++) {
                    for (int x=0; x<sizeX; x++) {
                        diffMap(y,x) = fabs(intMap7(y,x) - intMap3(y,x));
                    }
                }
            }
            else {
                for (int y=0; y<sizeY; y++) {
                    for (int x=0; x<sizeX; x++) {
                        diffMap(y,x) = fabs(intMap3(y,x));
                    }
                }
            }
            if (opmode & SaveMaps) {
                sprintf(fName, "%s_diff_%010d.gz", outfilename, i+1);
                if (diffMap.Write(fName))
                    cerr << "Error writing diff map " << fName << endl;
                else
                    cerr << fName << " written" << endl;
            }

            // sum of the pixel value and the 8-size neighbor's
            AgileMap sumMap;
            sumMap.ResizeTo(sizeY, sizeX);
            sumMap.SetElements(0.);
            for (int y=0; y<sizeY; y++) {
                for (int x=0; x<sizeX; x++) {
                    if ((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY) < radius*radius) {
                        sumMap(y,x) = diffMap(  y,   x) +
                                      diffMap(  y, x+1) +
                                      diffMap(y+1, x+1) +
                                      diffMap(y+1,   x) +
                                      diffMap(y+1, x-1) +
                                      diffMap(  y, x-1) +
                                      diffMap(y-1, x-1) +
                                      diffMap(y-1,   x) +
                                      diffMap(y-1, x+1);
                    }
                }
            }
            if (opmode & SaveMaps) {
                sprintf(fName, "%s_sumNeigh_%010d.gz", outfilename, i+1);
                if (sumMap.Write(fName))
                    cerr << "Error writing diff map " << fName << endl;
                else
                    cerr << fName << " written" << endl;
            }

            // get the max sum
            double max = 0.;
            int maxx = -1;
            int maxy = -1;
            for (int y=0; y<sizeY; y++) {
                for (int x=0; x<sizeX; x++) {
                    if(sumMap(y, x) > max) {
                        max = sumMap(y, x);
                        maxx = x;
                        maxy = y;
                    }
                }
            }
            os << i << " " << max << " " << maxx << " " << maxy << endl;
        }
    }
}
