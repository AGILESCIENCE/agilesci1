/*
 * Copyright (c) 2016
 *     Andrea Zoli (IASF-Bologna)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>

using std::cout;
using std::endl;

const char* startString = {
"################################################################\n"
"###               Task AG_select5 v1.0.0 - A.Z.              ###"
};

const char* endString = {
"### Task AG_select5 exiting ................................ ###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "indexfile", "Index file name" },
    { PilString, "timelist", "Time intervals file name" },
    { PilReal, "tmin", "Initial time (sec)" },
    { PilReal, "tmax", "Final time (sec)" },
    { PilNone, "", "" }
};

int main(int argc, char *argv[])
{
    cout << startString << endl;

    const char *expr = "";
    if (argc >= 7)
        expr = argv[6];

    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

    Intervals intervals;
    double tmin = params["tmin"];
    double tmax = params["tmax"];
    if (!eval::LoadTimeList(params["timelist"], intervals, tmin, tmax)) {
        cerr << "Error loading timelist file '" << params["timelist"].GetStr() << "'" << endl;
        return EXIT_FAILURE;
    }

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();

    cout << "Selecting the events.." << endl;
    const char *selectionFilename = params["outfile"];
    char *indexfile = (char*) params["indexfile"].GetStr();
    if (indexfile && indexfile[0]=='@')
        indexfile++;

    ostringstream oss;
    oss << selection::TimesExprString(intervals);
    if (strcmp(expr, "") != 0)
        oss << " && " << expr;
    cout << "Selection expr: " << oss.str() << endl;
    char templateFilename[FLEN_FILENAME];
    tmpnam(templateFilename);
    int status = selection::MakeSelection(indexfile, intervals, oss.str(), selectionFilename, templateFilename);
    FitsFile tindexfile(templateFilename);
    tindexfile.Delete();
    if (status==-118) {
        cout << endl << "AG_select5......................no matching events found" << endl;
        cout << endString << endl;
        return 0;
    }
    else if (status != 0) {
        cout << endl << "AG_select5......................selection failed" << endl;
        cout << endString << endl;
        return 0;
    }
    cout << endString << endl;

    return status;
}
