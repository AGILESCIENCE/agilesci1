/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

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
"### Task AG_expmapgen5 v1.3.0 - A.C., A.P.,  A.T., A.B., A.Z. ###"
};

const char* endString = {
"### Task AG_expmapgen5 exiting ............................... ###\n"
"##################################################################"
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
