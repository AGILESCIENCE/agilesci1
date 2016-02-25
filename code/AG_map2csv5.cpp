////////////////////////////////////////////////////////////////////////////////
// Date: Feb 2016
// Authors: Andrea Zoli (IASF_Bologna)
//
// INPUT
//       An Agile Map
//
// OUTPUT
//       A scv with a the list of "l-coord b-coord pixel-value" for each pixel of the Map
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <PilParams.h>
#include <AgileMap.h>

using std::cout;
using std::cerr;
using std::endl;

const char* startString = {
"###################################################\n"
"### Task AG_map2csv5 started ######################"
};

const char* endString = {
"### Task AG_map2csv5 ended successfully ###########\n"
"###################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "inmap", "Input map filename" },
    { PilString, "outfile", "Output filename" },
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
    const char *outFilename = params["outfile"];

    AgileMap map;
    if (map.Read(inMapFilename)) {
        cerr << "Error accessing file " << inMapFilename << endl;
        return EXIT_FAILURE;
    }

//    cout << "Center: " << map.GetMapCenterL() << " " << map.GetMapCenterB() << endl;
//    cout << "map.Dim(0) " << map.Dim(0) << " " << " map.Dim(1) " << map.Dim(1) << endl;
    std::ofstream ofs(outFilename);
    for (int y=0; y<map.Dim(0); ++y) {
        for (int x=0; x<map.Dim(1); ++x) {
            ofs << map.l(y, x) << " " << map.b(y, x) << " " << map(y, x) << endl;
        }
    }

    cout << endString << endl;
    return EXIT_SUCCESS;
}
