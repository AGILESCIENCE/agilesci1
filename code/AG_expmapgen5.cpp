////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_expmapgen
//       Authors: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
//       Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
//       Tomaso Contessi (Nuove Idee sas)
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
"################################################################\n"
"###   AG_expmapgen B25 v1.3.0 - A.C., T.C., A.T., A.B., A.Z. ###\n"
"################################################################\n"
};

const char* endString = {
"#################################################################\n"
"###  AG_expmapgen B25 exiting  .............................. ###\n"
"#################################################################\n"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "logfile", "Grid log index file name" },
    { PilString, "sarFileName", "Effective area file name" },
    { PilString, "edpFileName", "Energy dispersion file name" },
    { PilString, "maplist", "Maplist file name" },
    { PilString, "timelist", "Time intervals file name" },
    { PilReal, "mdim", "Size of Map (degrees)" },
    { PilReal, "mres", "Bin size (degrees)" },
    { PilReal, "la", "Longitude of map center (galactic)" },
    { PilReal, "ba", "Latitude of map center (galactic)" },
    { PilReal, "lonpole", "Rotation of map (degrees)" },
    { PilReal, "albrad", "Radius of earth albedo (degrees)" },
    { PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
    { PilReal, "roll_tol", "Roll tolerance (degrees)" },
    { PilReal, "earth_tol", "Earth tolerance (degrees)" },
    { PilInt, "phasecode", "Orbital phase code" },
    { PilString, "projection", "Projection (ARC or AIT)" },
    { PilReal, "binstep", "Bin step size" },
    { PilInt, "timestep", "LOG file step size" },
    { PilReal, "index", "Spectral index" },
    { PilReal, "tmin", "Initial time (sec)" },
    { PilReal, "tmax", "Final time (sec)" },
    { PilReal, "emin", "Minimum energy" },
    { PilReal, "emax", "Maximum energy" },
    { PilReal, "fovradmin", "Min radius of field of view (degrees)" },
    { PilReal, "fovradmax", "Max radius of field of view (degrees)" },
    { PilNone, "", "" }
};

int main(int argc, char *argv[])
{
    cout << startString << endl;

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

    cout << "INTERVALS N=" << intervals.Count() << ":" << endl;
    for (int i=0; i<intervals.Count(); ++i)
        cout << "   " << intervals[i].String() << endl;

    cout << "Selecting the events.." << endl;
    char selectionFilename[FLEN_FILENAME];
    char templateFilename[FLEN_FILENAME];
    tmpnam(selectionFilename);
    tmpnam(templateFilename);
    char *logfile = (char*) params["logfile"].GetStr();
    if (logfile && logfile[0]=='@')
        ++logfile;
    string logExpr = selection::LogExprString(intervals, params["phasecode"], params["timestep"]);
    int status = selection::MakeSelection(logfile, intervals, logExpr, selectionFilename, templateFilename);
    if (status != 0 && status != -118) {
        cout << endl << "AG_expmapgen5......................selection failed" << endl;
        cout << endString << endl;
        FitsFile sfile(selectionFilename);
        sfile.Delete();
        FitsFile tfile(templateFilename);
        tfile.Delete();
        return 0;
    }

    vector< vector<double> > exposures;
    status = eval::EvalExposure(params["outfile"], params["sarFileName"], params["edpFileName"],
                      params["maplist"], params["projection"], params["mdim"], params["mres"],
                      params["la"], params["ba"], params["lonpole"], params["albrad"],
                      params["y_tol"], params["roll_tol"], params["earth_tol"], params["phasecode"],
                      params["binstep"], params["timestep"], params["index"], params["tmin"], params["tmax"],
                      params["emin"], params["emax"], params["fovradmin"], params["fovradmax"],
                      selectionFilename, templateFilename, intervals, exposures, true);
    FitsFile sfile(selectionFilename);
    sfile.Delete();
    FitsFile tfile(templateFilename);
    tfile.Delete();

    if (status) {
        cout << "AG_expmapgen5..................... exiting AG_expmapgen ERROR:" << endl;
        cout << endString << endl;
        fits_report_error(stdout, status);
        return status;
    }
    cout << endString << endl;

    return status;
}
