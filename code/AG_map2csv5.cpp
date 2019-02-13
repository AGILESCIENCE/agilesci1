////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_map2csv
//		 Feb 2016
//       Author: Andrea Zoli (INAF/IASF Bologna)
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
//       Copyright (C) 2005-2019 AGILE Team. All rights reserved.
/*
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
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
"###    AG_map2csv  B25 v1.2.0 - A.Z.            ###\n"
"###################################################\n"
};

const char* endString = {
"###################################################\n"
"###      AG_map2csv5 ended successfully ###########\n"
"###################################################\n"
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
