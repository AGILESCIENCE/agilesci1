////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       AG_kernsim5
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

#include <TRandom3.h>
#include <RoiMulti5.h>
#include <PilParams.h>

using std::cout;
using std::cerr;
using std::endl;

const PilDescription c_params[] = {
    { PilInt,    "opmode", "Operation Mode" },
    { PilInt,    "nruns", "Number of runs" },
    { PilInt,    "seed", "Seed" },
    { PilInt,    "radius", "Radious of the kernel, radius 1 for a 3x3 matrix" },
    { PilString, "sarfile", "SAR file name" },
    { PilString, "edpfile", "EDP file name" },
    { PilString, "psdfile", "PSD file name" },
    { PilString, "maplistsim", "Map list for simulation" },
    { PilString, "srclistsim", "Source list for simulation" },
    { PilString, "outfile", "Output file prefix" },
    { PilNone,   "", "" }
};

class KernSimParams: public PilParams {
public:
    KernSimParams() : PilParams(c_params) {}
    void Print() { PilParams::Print(c_params); }
};


enum { Analysis=1, SaveMaps=2 };

class AppScreen {
public:
    AppScreen() {
        cout << "#################################################################"<< endl;
        cout << "######## AG_kernsim5  v.1.0 - 11/01/2016 - A.Z., A.B. ############"<< endl;
        cout << "#################################################################"<< endl;
        cout << "#################################################################"<< endl;
        cout << endl << "RoiMulti version: " << ROIMULTI_VERSION << " - " << ROIMULTI_DATE << endl;
    }

    ~AppScreen() {
        cout << "#################################################################"<< endl;
        cout << "##########  Task AG_kernsim5............ exiting #################"<< endl;
        cout << "#################################################################"<< endl;
    }
};

int main(int argc,char **argv) {
    AppScreen appScreen;

    KernSimParams mPars;
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
    const char* maplistsimname = mPars["maplistsim"];
    const char* srclistsim = mPars["srclistsim"];
    const char* outfilename = mPars["outfile"];

    if (seed)
        SetSeed(seed);

    // compute simulated intensity maps
    MapList maplistsim;
    if (!maplistsim.Read(maplistsimname))
        return -1;

    MapData mapData;
    if (!mapData.Load(maplistsim, true)) {
        cerr << "Error reading the simulation map list " << endl;
        return -1;
    }

    SourceDataArray srcSimArr = ReadSourceFile(srclistsim);

    char fName[256];
    int sizeX = 0;
    int sizeY = 0;
    AgileMap sumMap;
    for (int i=0; i<nruns; ++i) {
        RoiMulti roiMulti;
        if (!roiMulti.SetPsf(psdfilename, sarfilename, edpfilename)) {
            cerr << "ERROR setting PSF data" << endl;
            return -1;
        }

        cout << endl << "AG_kernsim5 loop #" << i+1 << endl << endl;

        mapData.MapCoeff::Load(maplistsim);
        roiMulti.SetMaps(mapData);
        AgileMap* simArr = roiMulti.NewSimulationArray(srcSimArr);

        if (opmode & SaveMaps) {
            sprintf(fName, "%s_%010d.cts.gz", outfilename, i+1);
            if (simArr[0].Write(fName))
                cerr << "Error writing simulated counts map." << fName << endl;
            else
                cerr << fName << " written" << endl;
        }
        if (opmode & Analysis) {
            const AgileMap& ctsMap = simArr[0];
            const AgileMap& expMap = mapData.ExpMap(0);

            // first loop analysis init
            if (i == 0) {
                sizeY = expMap.Dim(0);
                sizeX = expMap.Dim(1);
            }

            // compute intensity maps
            AgileMap intMap;
            intMap.ResizeTo(sizeY, sizeX);
            sumMap.ResizeTo(sizeY, sizeX);
            for (int y=0; y<sizeY; y++) {
                for (int x=0; x<sizeX; x++) {
                    double intensity;
                    if(expMap(y,x) != 0)
                        intensity = ctsMap(y,x) / expMap(y,x);
                    else
                        intensity = 0.0;
/*                    if(intensity < 0 || intensity > 0.012)
                        intensity = 0.0;*/
                    intMap(y,x) = intensity;
                    sumMap(y,x) += intensity;
                }
            }

            if (opmode & SaveMaps) {
                sprintf(fName, "%s_%010d.int.gz", outfilename, i+1);
                if (intMap.Write(fName))
                    cerr << "Error writing simulated intensity map." << fName << endl;
                else
                    cerr << fName << " written" << endl;
            }
        }
    }

    // compute the kern intensity map
    // save the maximum to file
    // compute the image sum/integral
    double max = 0.;
    int maxx = -1;
    int maxy = -1;
    AgileMap kernMap;
    kernMap.ResizeTo(sizeY, sizeX);
    for (int y=0; y<sizeY; y++) {
        for (int x=0; x<sizeX; x++) {
            kernMap(y, x) = sumMap(y,x) / nruns;
            if(fabs(kernMap(y, x)) > fabs(max)) {
                max = kernMap(y, x);
                maxx = x;
                maxy = y;
            }
        }
    }
    sprintf(fName, "%s_kern.gz", outfilename);
    if (kernMap.Write(fName))
        cerr << "Error writing kern map." << fName << endl;
    else
        cerr << fName << " written" << endl;

    // save the kernel
    sprintf(fName, "%s_kern.csv", outfilename);
    std::ofstream os(fName);
    double sum = 0.0;
    for (int x=maxx-radius; x<=maxx+radius; x++) {
        for (int y=maxy-radius; y<=maxy+radius; y++) {
            os << kernMap(y, x);
            if (y < maxy + radius)
                os << ", ";
            sum += kernMap(y, x);
        }
        os << endl;
    }
    os.close();

    // save the kernel normalized (integral = 1)
    sprintf(fName, "%s_kern_norm.csv", outfilename);
    os.open(fName);
    for (int x=maxx-radius; x<=maxx+radius; x++) {
        for (int y=maxy-radius; y<=maxy+radius; y++) {
            os << kernMap(y, x) / sum;
            if (y < maxy + radius)
                os << ", ";
        }
        os << endl;
    }
    os.close();

    return EXIT_SUCCESS;
}
