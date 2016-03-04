////////////////////////////////////////////////////////////////////////////////
// Date: Jan 2016
// Authors: Andrea Zoli (IASF_Bologna)
//
// INPUT
//       A map, center of a circle (l,b) in degrees (galactic) and its radius.
//
// OUTPUT
//       The input map with a circle drawn. Every circle value is a 1.0.
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <PilParams.h>
#include <AgileMap.h>

using std::cout;
using std::cerr;
using std::endl;

const char* startString = {
"###################################################\n"
"### Task AG_circle5 v1.1.0 - A.Z.               ###"
};

const char* endString = {
"### Task AG_circle5 ended successfully ############\n"
"###################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "inmap", "Input map filename" },
    { PilString, "outmap", "Output map filename" },
    { PilReal,   "l", "Circle center longitude l in degrees (galactic)" },
    { PilReal,   "b", "Circle center latitude b in degrees (galactic)" },
    { PilReal,   "radius", "Radius of the circle in degrees" },
    { PilNone,   "",   "" }
};

int main(int argc, char *argv[]) {
    cout << startString << endl;

    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();

    const char *inMapFilename = params["inmap"];
    const char *outMapFilename = params["outmap"];

    AgileMap map;
    if (map.Read(inMapFilename)) {
        cerr << "Error accessing file " << inMapFilename << endl;
        return EXIT_FAILURE;
    }

    double l = params["l"];
    double b = params["b"];
    double radius = params["radius"];
    int row, col;
    map.GetRowCol(l, b, &row, &col);
    cout << "Circle coordinates:" << endl;
    cout << "l " << l   << "   b " << b << endl;
    cout << "y " << row << "   x " << col << endl;
    for (int y=0; y<map.Dim(0); ++y) {
        for (int x=0; x<map.Dim(1); ++x) {
            double dist = map.SrcDist(y, x, l, b);
            if(dist < radius)
                map(y, x) = 1.;
            else
                map(y, x) = 0.;
        }
    }

    if(map.Write(outMapFilename))
        return EXIT_FAILURE;

    cout << endString << endl;
    return EXIT_SUCCESS;
}
