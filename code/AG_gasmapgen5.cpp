////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_gasmapgen
//       Author: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
//       Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
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
#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>

using std::cout;
using std::endl;
using std::vector;

const char* startString = {
"#################################################################\n"
"###    AG_expmapgen B25 v1.3.0 - A.C. A.P. A.T. A.B. A.Z.     ###\n"
"#################################################################\n"
};

const char* endString = {
"##################################################################\n"
"###   AG_expmapgen B25 exiting ............................... ###\n"
"##################################################################\n"
};

const PilDescription paramsDescr[] = {
    { PilString, "expfile", "Exposure filename" },
    { PilString, "outfile", "Output filename" },
    { PilString, "diffusefile", "Diffuse model filename" },
    { PilString, "hiresdiffusefile", "High res diffuse model filename" },
    { PilNone, "", "" }
};

int main(int argc, char *argv[])
{
    cout << startString << endl;

    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();

    int status = eval::EvalGas(params["outfile"], params["expfile"], params["diffusefile"], params["hiresdiffusefile"]);
    if (status) {
        cout << "AG_gasmapgen5..................... exiting AG_gasmapgen ERROR:" << endl;
        cout << endString << endl;
        fits_report_error(stdout, status);
        return status;
    }
    cout << endString << endl;

    return status;
}
