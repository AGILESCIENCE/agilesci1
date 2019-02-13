////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG summapgen
//       Feb 2016
//		 Author: Andrea Zoli (INAF/IASF Bologna)
//
// INPUT
//       A maplist4, with multiple cts and exposure maps.
//
// OUTPUT
//       Gives a cts and an exp map computed respectively as the sum of cts maps and the sum of exp maps.
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
#include <cstdlib>
#include <PilParams.h>
#include <AgileMap.h>
#include <AlikeData5.h>

using std::cout;
using std::cerr;
using std::endl;

const char* startString = {
"###################################################\n"
"###    AG_summapgen B25 v1.2.0 - A.Z.           ###\n"
"###################################################\n"
};

const char* endString = {
"###################################################\n"
"###   AG_summapgen B25 ended successfully #########\n"
"###################################################\n"
};

const PilDescription paramsDescr[] = {
    { PilString, "maplist", "Input maplist" },
    { PilString, "outprefix", "Output name prefix for the maps" },
    { PilString, "operationmode", "Operation mode: sum (it adds the maps), sub (it subtracts from the first map the others map)" },
    { PilNone, "", "" }
};

int main(int argc, char *argv[]) {
    cout << startString << endl;

    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();

    MapList maplist;
    int mapCount = maplist.Read(params["maplist"]);
    if (!mapCount) {
        cerr << "File " << params.GetStrValue("maplist") << " missing or empty." << endl;
        return EXIT_FAILURE;
    }
    MapData mapData;
    if (!mapData.Load(maplist)) {
        cerr << "Error loading map data from " << params.GetStrValue("maplist") << "." << endl;
        return EXIT_FAILURE;
    }

    std::string operationMode = std::string(params["operationmode"]);
    short int om;
    if(operationMode == "sum")
        om = +1;
    else if(operationMode == "sub")
        om = -1;
    else
    {
        cerr << "Operation mode not correct. Possible values: [sum, sub]." << endl;
        return EXIT_FAILURE;
    }


    AgileMap sumCts = mapData.CtsMap(0);
    AgileMap sumExp = mapData.ExpMap(0);
    double emin = sumCts.GetEmin();
    double emax = sumCts.GetEmax();
    double fovmin = sumCts.GetFovMin();
    double fovmax = sumCts.GetFovMax();
    double tstart = sumCts.GetTstart();
    double tstop = sumCts.GetTstop();
    for(int i=1; i<mapData.Count(); i++) {
        AgileMap otherCts = mapData.CtsMap(i);
        for (int y=0; y<sumCts.Dim(0); ++y)
            for (int x=0; x<sumCts.Dim(1); ++x)
            {
                sumCts(y, x) += (om) * otherCts(y, x);
                if(sumCts(y, x) < 0)
                   sumCts(y, x) = 0;
            }

        AgileMap otherExp = mapData.ExpMap(i);
        for (int y=0; y<sumExp.Dim(0); ++y)
            for (int x=0; x<sumExp.Dim(1); ++x)
            {
                sumExp(y, x) += (om) * otherExp(y, x);
                if (sumExp(y, x) < 0 )
                    sumExp(y, x) = 0;
            }
        if(otherCts.GetEmin() < emin)
            emin = otherCts.GetEmin();
        if(otherCts.GetEmax() > emax)
            emax = otherCts.GetEmax();
        if(otherCts.GetFovMin() < fovmin)
            fovmin = otherCts.GetFovMin();
        if(otherCts.GetFovMax() > fovmax)
            fovmax = otherCts.GetFovMax();
        if(otherCts.GetTstart() < tstart)
            tstart = otherCts.GetTstart();
        if(otherCts.GetTstop() > tstop)
            tstop = otherCts.GetTstop();
    }
    sumCts.SetEnergy(emin, emax);
    sumExp.SetEnergy(emin, emax);
    sumCts.SetFov(fovmin, fovmax);
    sumExp.SetFov(fovmin, fovmax);
    sumCts.SetTT(tstart, tstop);
    sumExp.SetTT(tstart, tstop);

    std::string outCts = std::string(params["outprefix"])+".cts.gz";
    if(sumCts.WriteWithAllMetadata(outCts.c_str()))
        return EXIT_FAILURE;

    std::string outExp = std::string(params["outprefix"])+".exp.gz";
    if(sumExp.WriteWithAllMetadata(outExp.c_str()))
        return EXIT_FAILURE;

    cout << endString << endl;
    return EXIT_SUCCESS;
}
