/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>

#include "EVTFilter.h"
#include "LOGFilter.h"

using std::cout;
using std::endl;

const char* startString = {
"################################################################\n"
"###                   Task AG_ap5 v1.0.0 - A.B. A.Z.         ###"
};

const char* endString = {
"### Task AG_ap5 exiting .................................... ###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "logfile", "Grid log index file name" },
    { PilString, "evtfile", "Event file index file name" },
    { PilString, "sarFileName", "Effective area file name" },
    { PilString, "edpFileName", "Energy dispersion file name" },
    { PilString, "timelist", "Time intervals file name" },
    { PilReal, "mres", "Bin size (degrees)" },
    { PilReal, "la", "Longitude of map center (galactic)" },
    { PilReal, "ba", "Latitude of map center (galactic)" },
    { PilReal, "lonpole", "Rotation of map (degrees)" },
    { PilReal, "albrad", "Radius of earth albedo (degrees)" },
    { PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
    { PilReal, "roll_tol", "Roll tolerance (degrees)" },
    { PilReal, "earth_tol", "Earth tolerance (degrees)" },
    { PilInt, "phasecode", "Orbital phase code" },
    { PilInt, "timestep", "LOG file step size" },
    { PilReal, "index", "Spectral index" },
    { PilReal, "tmin", "Initial time (sec)" },
    { PilReal, "tmax", "Final time (sec)" },
    { PilReal, "emin", "Minimum energy" },
    { PilReal, "emax", "Maximum energy" },
    { PilReal, "fovradmin", "Min radius of field of view (degrees)" },
    { PilReal, "fovradmax", "Max radius of field of view (degrees)" },
    { PilInt, "filtercode", "Event filter code" },
    { PilReal, "timeslot", "Time slot" },
    { PilInt, "usetelem", "Use telemetry" },
    { PilNone, "", "" }
};

int main(int argc, char *argv[])
{
    cout << startString << endl;

    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

    int timestep = params["timestep"];
    Intervals intervals;
    double tmin = params["tmin"];
    double tmax = params["tmax"];
    if (!eval::LoadTimeList(params["timelist"], intervals, tmin, tmax)) {
        cerr << "Error loading timelist file '" << params["timelist"].GetStr() << "'" << endl;
        return EXIT_FAILURE;
    }

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();
    double mdim = params["mres"];
    mdim = mdim * sqrt(2);
    double radius = params["mres"];
    double binstep = 1.0;
    const char *projection = "ARC";
    cout << "### Radius for evt: " << radius << " - mdim for exp: " << mdim << endl;
    cout << "### Binstep: " << binstep << endl;
    cout << "### Projection: " << projection << endl;

    //area calculation
    double pixel1 = DEG2RAD * DEG2RAD * fabs(mdim * mdim);
    double areapixel =  pixel1 * Sinaa(DEG2RAD*45.);

    cout << "### Area pixel " << areapixel << endl;

    cout << "INTERVALS N=" << intervals.Count() << ":" << endl;
    for (int i=0; i<intervals.Count(); i++)
        cout << "   " << intervals[i].String() << endl;

    int usetelem = params["usetelem"];
    EVTFilter* evtfilter = 0;
    LOGFilter* logfilter = 0;
    char selectionLogFilename[FLEN_FILENAME];
    char templateLogFilename[FLEN_FILENAME];
    char selectionEvtFilename[FLEN_FILENAME];
    char templateEvtFilename[FLEN_FILENAME];

    char *evtfile = (char*) params["evtfile"].GetStr();
    if (evtfile && evtfile[0]=='@')
        evtfile++;

    char *logfile = (char*) params["logfile"].GetStr();
    if (logfile && logfile[0]=='@')
        logfile++;

    int status;

    if(usetelem == 0) {
      cout << "Selecting the events.." << endl;

      tmpnam(selectionLogFilename);
      tmpnam(templateLogFilename);

      string logExpr = selection::LogExprString(intervals, params["phasecode"], params["timestep"]);
      cout << logExpr << endl;
      status = selection::MakeSelection(logfile, intervals, logExpr, selectionLogFilename, templateLogFilename);
      if (status==-118) {
          cout << endl << "AG_ap5......................no matching events found" << endl;
          cout << endString << endl;
          return 0;
      }
      else if (status != 0) {
          cout << endl << "AG_ap5......................selection failed" << endl;
          cout << endString << endl;
          return 0;
      }

      cout << "Selecting the events.." << endl;

      tmpnam(selectionEvtFilename);
      tmpnam(templateEvtFilename);

      string evtExpr = selection::EvtExprString(intervals, params["emin"], params["emax"],
                                      params["albrad"], params["fovradmax"], params["fovradmin"],
                                      params["phasecode"], params["filtercode"]);
      status = selection::MakeSelection(evtfile, intervals, evtExpr, selectionEvtFilename, templateEvtFilename);
      if (status==-118) {
          cout << endl << "AG_ap5......................no matching events found" << endl;
          cout << endString << endl;
          return 0;
      }
      else if (status != 0) {
          cout << endl << "AG_ap5......................selection failed" << endl;
          cout << endString << endl;
          return 0;
      }
    } else {
      evtfilter = new EVTFilter(evtfile);
      logfilter = new LOGFilter(logfile, timestep);
    }

    const char *outfile = params["outfile"];
    std::ofstream expText(outfile);
    expText.setf(ios::fixed);
    double exp;
    double beginTime = tmin;
    double deltaT = params["timeslot"];
    double endTime = beginTime+deltaT;
    if (endTime > tmax)
        endTime = tmax;
    cout.setf(ios::fixed);
    cout << std::setprecision(2);
    cout << "***** " << beginTime << " " << endTime << " " << deltaT << endl << endl;
    double totalExp = 0;
    int totalCounts = 0;
    Interval timeSlot;
    do {
#ifdef DEBUG
        cout << "Time slot beginTime: " << beginTime << " endTime: " << endTime << endl;
#endif
        timeSlot.Set(beginTime, endTime);
        Intervals intervalSlots = Intersection(intervals, timeSlot);
        if (intervalSlots.Count()) {
            cout << "Selected slots:" << endl;
            for (int i=0; i<intervalSlots.Count(); i++)
                cout << "   " << intervalSlots[i].Start() << " " << intervalSlots[i].Stop() << endl;

            vector< vector<double> > exposures;
            if(usetelem == 0) {
              status = eval::EvalExposure("None", params["sarFileName"], params["edpFileName"],
                               "None", projection, mdim, mdim, params["la"], params["ba"],
                               params["lonpole"], params["albrad"], params["y_tol"], params["roll_tol"],
                               params["earth_tol"], params["phasecode"], binstep, params["timestep"],
                               params["index"], tmin, tmax, params["emin"],
                               params["emax"], params["fovradmin"], params["fovradmax"],
                               selectionLogFilename, templateLogFilename, intervalSlots, exposures, false, 0);
            } else {
              int phasecode = params["phasecode"];
              logfilter->query( tmin, tmax, phasecode );
              status = eval::EvalExposure("None", params["sarFileName"], params["edpFileName"],
                             "None", projection, mdim, mdim, params["la"], params["ba"],
                             params["lonpole"], params["albrad"], params["y_tol"], params["roll_tol"],
                             params["earth_tol"], params["phasecode"], binstep, params["timestep"],
                             params["index"], tmin, tmax, params["emin"],
                             params["emax"], params["fovradmin"], params["fovradmax"],
                             0, 0, intervalSlots, exposures, false, logfilter);
            }
			//TBW
			/*
			vector<double>  dist_pl_earth;
			vector<double>  dist_pl_source;
            status = eval::GetPLDirection("None", params["sarFileName"], params["edpFileName"],
                               "None", projection, mdim, mdim, params["la"], params["ba"],
                               params["lonpole"], params["albrad"], params["y_tol"], params["roll_tol"],
                               params["earth_tol"], params["phasecode"], binstep, params["timestep"],
                               params["index"], tmin, tmax, params["emin"],
                               params["emax"], params["fovradmin"], params["fovradmax"],
                               selectionLogFilename, templateLogFilename, intervalSlots, dist_pl_earth, dist_pl_source);
            */
      			vector<int>  counts;
            if(usetelem == 0) {
      			  status = eval::EvalCountsInRadius("None", tmin, tmax, radius,
      						   params["la"], params["ba"], params["lonpole"],
      						   params["emin"], params["emax"], params["fovradmax"],
      						   params["fovradmin"], params["albrad"], params["phasecode"],
      						   params["filtercode"], selectionEvtFilename, templateEvtFilename,
      						   intervalSlots, counts, 0);
            } else {
              double mdim = params["mdim"];
            	double la = params["la"];
            	double ba = params["ba"];
            	//evtfilter->setPostfilter1(mdim, la, ba);

              int phasecode = params["phasecode"];
            	int filtercode =  params["filtercode"];
            	double emin = params["emin"];
            	double emax = params["emax"];
            	double albrad =  params["albrad"];
            	double fovradmin =  params["fovradmin"];
            	double fovradmax =  params["fovradmax"];
              evtfilter->query(tmin, tmax, phasecode, filtercode, emin, emax, albrad, fovradmin, fovradmax );

              status = eval::EvalCountsInRadius("None", tmin, tmax, radius,
                    params["la"], params["ba"], params["lonpole"],
                    params["emin"], params["emax"], params["fovradmax"],
                    params["fovradmin"], params["albrad"], params["phasecode"],
                    params["filtercode"], 0, 0,
                    intervalSlots, counts, evtfilter);

            }
            double slotExp = 0;
            int slotCounts = 0;
            for (int slot=0; slot<intervalSlots.Count(); slot++) {
                slotExp += exposures[slot][0]; // the map is 1x1
                slotCounts += counts[slot];
            }
			//slotExp in cm2 s sr
			//output in cm2 s
            if(status == 0) {
                expText << std::setprecision(1);
                expText << beginTime << " " << endTime << " ";
                expText << std::setprecision(2);
                expText << slotExp / areapixel  << " " << slotCounts << " ";
                /*
                for (int slot=0; slot<intervalSlots.Count(); slot++) {
                	expText << dist_pl_earth[slot] << " " << dist_pl_source[slot] << " ";
                }
                */
                totalExp += slotExp;
                totalCounts += slotCounts;
                expText << endl;
            }
            else if(status == -118)
                break;
        }
        else
            cout << "No intervals selected" << endl;

        beginTime = endTime;
        endTime += deltaT;
        if (tmax < endTime)
            endTime = tmax;
    } while (beginTime < tmax);
    expText.close();
    cout << "Total Counts: " << totalCounts << endl;
    cout << "Total Exposure [cm2 s sr]: " << totalExp << endl;

    if(usetelem == 0) {
      FitsFile slogfile(selectionLogFilename);
      slogfile.Delete();
      FitsFile tlogfile(templateLogFilename);
      tlogfile.Delete();
      FitsFile sevtfile(selectionEvtFilename);
      sevtfile.Delete();
      FitsFile tevtfile(templateEvtFilename);
      tevtfile.Delete();
    }
    if (status == -118) {
        cout << endl << "AG_ap5......................no matching events found" << endl;
    }
    else if (status != 0) {
        cout << endl << "AG_ap5...................... exiting with ERROR:"<< endl;
        if(usetelem == 0)
          fits_report_error(stdout, status);
    }
    cout << endString << endl;

    return status;
}
