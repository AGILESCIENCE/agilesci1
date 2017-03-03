/*
 * Copyright (c) 2017
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
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

using std::cout;
using std::endl;

const char* startString = {
"################################################################\n"
"###                   Task AG_lm5 v1.0.0 - A.B.              ###"
};

const char* endString = {
"### Task AG_lm5 exiting .................................... ###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "logfile", "Grid log index file name" },
    { PilString, "evtfile", "Event file index file name" },
    { PilString, "sarFileName", "Effective area file name" },
    { PilString, "edpFileName", "Energy dispersion file name" },
    { PilString, "timelist", "Time intervals file name" },
    { PilReal, "lonpole", "Rotation of map (degrees)" },
    { PilReal, "albrad", "Radius of earth albedo (degrees)" },
    { PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
    { PilReal, "roll_tol", "Roll tolerance (degrees)" },
    { PilReal, "earth_tol", "Earth tolerance (degrees)" },
    { PilInt, "phasecode", "Orbital phase code" },
    { PilInt, "timestep", "LOG file step size" },
    { PilReal, "index", "Spectral index" },
    { PilReal, "emin", "Minimum energy" },
    { PilReal, "emax", "Maximum energy" },
    { PilReal, "fovradmin", "Min radius of field of view (degrees)" },
    { PilReal, "fovradmax", "Max radius of field of view (degrees)" },
    { PilInt, "filtercode", "Event filter code" },
    { PilReal, "t0", "T0 (TT)" },
    { PilReal, "la", "Longitude of GRB centroid (galactic)" },
    { PilReal, "ba", "Latitude of GRB centroid (galactic)" },
    { PilReal, "radius", "LM radius of analysis" },   
    { PilReal, "t1s", "t1s (sec) - analysis of signal in [T0-t1s, T0+t2s]" },
    { PilReal, "t2s", "t2s (sec) - analysis of signal in [T0-t1s, T0+t2s]" },
    { PilReal, "t1b", "t1b (sec)" },
    { PilReal, "shiftt1b", "shift t1b (sec) - analysis of background 1 in [T0-t1s-shiftt1b-t1b, T0-t1s-shiftt1b[" },
    { PilReal, "t2b", "t2b (sec)" },
    { PilReal, "shiftt2b", "shift t2b (sec) - analysis of background 2 in ]T0+t2s+shiftt2b, T0+t2s+shiftt2b+t2b]" },
    { PilNone, "", "" }
};


int EvalExpAndCounts(PilParams &params, double tmin, double tmax, int &countscalc, double &expcalc)
{
	
	int timestep = 1;
	
    Intervals intervals;
    if (!eval::LoadTimeList(params["timelist"], intervals, tmin, tmax)) {
        cerr << "Error loading timelist file '" << params["timelist"].GetStr() << "'" << endl;
        return EXIT_FAILURE;
    }

    cout << endl << "INPUT PARAMETERS:" << endl;
    params.Print();
    double radius = params["radius"]; //mres
    double mdim = radius;
    double binstep = 1.0;
    const char *projection = "ARC";
    cout << endl << "Mdim: " << mdim << endl;
    cout << "Binstep: " << binstep << endl;
    cout << "Projection: " << projection << endl;

    cout << "INTERVALS N=" << intervals.Count() << ":" << endl;
    for (int i=0; i<intervals.Count(); i++)
        cout << "   " << intervals[i].String() << endl;

    cout << "Selecting the events.." << endl;
    char selectionLogFilename[FLEN_FILENAME];
    char templateLogFilename[FLEN_FILENAME];
    tmpnam(selectionLogFilename);
    tmpnam(templateLogFilename);
    char *logfile = (char*) params["logfile"].GetStr();
    if (logfile && logfile[0]=='@')
        logfile++;
    string logExpr = selection::LogExprString(intervals, params["phasecode"], timestep);
    int status = selection::MakeSelection(logfile, intervals, logExpr, selectionLogFilename, templateLogFilename);
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
    char selectionEvtFilename[FLEN_FILENAME];
    char templateEvtFilename[FLEN_FILENAME];
    tmpnam(selectionEvtFilename);
    tmpnam(templateEvtFilename);
    char *evtfile = (char*) params["evtfile"].GetStr();
    if (evtfile && evtfile[0]=='@')
        evtfile++;
    string evtExpr = selection::EvtExprString(intervals, params["emin"], params["emax"],
                                    params["albrad"], params["fovradmax"], params["fovradmin"],
                                    params["phasecode"], params["filtercode"]);
    status = selection::MakeSelection(evtfile, intervals, evtExpr, selectionEvtFilename, templateEvtFilename);
    if (status==-118) {
        cout << endl << "AG_lm5......................no matching events found" << endl;
        cout << endString << endl;
        return 0;
    }
    else if (status != 0) {
        cout << endl << "AG_lm5......................selection failed" << endl;
        cout << endString << endl;
        return 0;
    }

    double beginTime = tmin;
   	double endTime = tmax;
    cout.setf(ios::fixed);
    cout << std::setprecision(2);
    cout << "***** " << beginTime << " " << endTime << " " <<  endl << endl;
    Interval timeSlot;
    timeSlot.Set(beginTime, endTime);
    
#ifdef DEBUG
    cout << "Time slot beginTime: " << beginTime << " endTime: " << endTime << endl;
#endif
	
	Intervals intervalSlots = Intersection(intervals, timeSlot);
	if (intervalSlots.Count()) {
		cout << "Selected slots:" << endl;
		for (int i=0; i<intervalSlots.Count(); i++)
			cout << "   " << intervalSlots[i].Start() << " " << intervalSlots[i].Stop() << endl;

		vector< vector<double> > exposures;
		status = eval::EvalExposure("None", params["sarFileName"], params["edpFileName"],
						   "None", projection, mdim, params["radius"], params["la"], params["ba"],
						   params["lonpole"], params["albrad"], params["y_tol"], params["roll_tol"],
						   params["earth_tol"], params["phasecode"], binstep, params["timestep"],
						   params["index"], tmin, tmax, params["emin"],
						   params["emax"], params["fovradmin"], params["fovradmax"],
						   selectionLogFilename, templateLogFilename, intervalSlots, exposures, false);

		vector<int>  counts;
		status = eval::EvalCountsInRadius("None", tmin, tmax, radius, 
						   params["la"], params["ba"], params["lonpole"],
						   params["emin"], params["emax"], params["fovradmax"],
						   params["fovradmin"], params["albrad"], params["phasecode"],
						   params["filtercode"], selectionEvtFilename, templateEvtFilename,
						   intervalSlots, counts);
			             
		expcalc = 0;
		countscalc = 0;
		for (int slot=0; slot<intervalSlots.Count(); slot++) {
			expcalc += exposures[slot][0]; // the map is 1x1
			countscalc += counts[slot]; 
		}

		
	}
	else
		cout << "No intervals selected" << endl;

    FitsFile slogfile(selectionLogFilename);
    slogfile.Delete();
    FitsFile tlogfile(templateLogFilename);
    tlogfile.Delete();
    FitsFile sevtfile(selectionEvtFilename);
    sevtfile.Delete();
    FitsFile tevtfile(templateEvtFilename);
    tevtfile.Delete();

    if (status == -118) {
        cout << endl << "AG_lm5......................no matching events found" << endl;
    }
    else if (status != 0) {
        cout << endl << "AG_lm5...................... exiting with ERROR:"<< endl;
        fits_report_error(stdout, status);
    }
    cout << endString << endl;

    return status;
}

int main(int argc, char *argv[])
{
    cout << startString << endl;
	int status = 0;
	
    PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;

	double t0 = params["t0"];
	double t1s = params["t1s"];
	double t2s = params["t2s"];
	double t1b = params["t1b"];
	double t2b = params["t2b"];
	double shiftt1b = params["shiftt1b"];
	double shiftt2b = params["shiftt2b"];
	
	double tmin = 0.0;
	double tmax = 0.0;
	int counts_s = 0;
	double exp_s = 0.0;
	int counts_b1 = 0;
	double exp_b1 = 0.0;
	int counts_b2 = 0;
	double exp_b2 = 0.0;
		
	const char *outfile = params["outfile"];
    std::ofstream resText(outfile);
    resText.setf(ios::fixed);
    
    /*
    double deltaT = params["timeslot"];
    double endTime = beginTime+deltaT;
    if (endTime > tmax)
        endTime = tmax;
	*/
	
	tmin = t0-t1s;
    tmax = t0+t2s;
    
    status = EvalExpAndCounts(params, tmin, tmax, counts_s, exp_s);
    
    if(status == 0) {
		resText << std::setprecision(1);
		resText << tmin << " " << tmax << " ";
		resText << std::setprecision(2);
		resText << counts_s << " " << exp_s << endl;
	}
	else if(status == -118)
		return status;

	tmin = t0-t1s-shiftt1b-t1b;
    tmax = t0-t1s-shiftt1b;
    
    status = EvalExpAndCounts(params, tmin, tmax, counts_b1, exp_b1);
    
    if(status == 0) {
		resText << std::setprecision(1);
		resText << tmin << " " << tmax << " ";
		resText << std::setprecision(2);
		resText << counts_b1 << " " << exp_b1 << endl;
	}
	else if(status == -118)
		return status;
    
    tmin = t0+t2s+shiftt2b;
    tmax = t0+t2s+shiftt2b+t2b;
    
    status = EvalExpAndCounts(params, tmin, tmax, counts_b2, exp_b2);
    
    if(status == 0) {
		resText << std::setprecision(1);
		resText << tmin << " " << tmax << " ";
		resText << std::setprecision(2);
		resText << counts_b2 << " " << exp_b2 << endl;
	}
	else if(status == -118)
		return status;
    
    int bkg = counts_b1 + counts_b2;
    int source = counts_s;
    int N_on = source + bkg;
    int N_off = bkg;
    double alpha = exp_s / (exp_b1 + exp_b2);
    resText << alpha << endl;
    double alp1 = alpha / (1 + alpha);
    double alp2 = alpha + 1;
    
    cout << "bkg rate " << bkg / (double) (exp_b1 + exp_b2) << endl;
    cout << "sig rate " << source / (double) exp_s << endl;
    
    if ((source > 0) and (bkg > 0)) {
    	double source1 = source;
    	double bkg1 = bkg;
        double L1 = pow(((source1 + bkg1) / source1) * alp1, source);
        double L2 = pow(((bkg1 + source1) / bkg1) / alp2, bkg);
        double L = L1 * L2;
        double S = sqrt(-2. * log(L));
        cout <<  "Li&Ma sigma " << S << endl;
    } else {
        cout << "Alpha: 0" << endl;
       	cout << "Li&Ma sigma 0" << endl;
    }
    /*
    beginTime = endTime;
	endTime += deltaT;
	if (tmax < endTime)
		endTime = tmax;
    */
    
    
    resText.close();
    
    return 0;
}

