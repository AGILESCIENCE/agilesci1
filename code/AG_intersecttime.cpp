/***************************************************************************
                          main.cpp
                             -------------------
    copyright            : (C) 2013 AGILE Team
    email                : bulgarelli@iasfbo.inaf.it
    contributors		 : Andrea Bulgarelli (IASF-Bologna)

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

//#include <LSSpectralAlgorithmNR.h>
//#include <TimedLoader.h>
//#include <SpectralAnalysis.h>
//#include <common_string.h>
#include <Intervals.h>

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//mainINTERSECTTIME
int main(int argc,char **argv) {

	if(argc <= 2) {
		cout << "0: time list of the good observations" << endl;
		cout << "1: time list to intersect" << endl;
		cout << "Timelist is a list of tstart tstop in TT" << endl;
		exit(0);
	}

    const char* intFileName = argv[1];
	
	Intervals intvs;
    
    intvs = ReadIntervalsNoUnion(intFileName);
    

    if (intvs.Count()>1) {
    	//cout << intvs.String() << endl;
        //cout << intvs.Count() << " intervals:" << endl;
        //for (int i=0; i<intvs.Count(); ++i)
        //    cout << "   " << String(intvs[i]) << endl;
    }
    
    const char* intFileName2 = argv[2];
    Intervals intvs2 = ReadIntervalsNoUnion(intFileName2);
    

	for(int i=0; i< intvs2.Count(); i++) {
		//cout << intervalSlots.String() << endl;
		//cout << "--" << endl;
		Interval inter = intvs2[i];
		Intervals intervalSlots = Intersection(intvs, inter);
		//cout << "-" << endl;
		if(intervalSlots.Count() > 0)
			cout << intervalSlots.String() << endl;
	}
	
	//cout << intervalSlots.Count() << " intervals intersection:" << endl;
    //for (int i=0; i<intervalSlots.Count(); ++i)
    //    cout << "   " << String(intervalSlots[i]) << endl;
    //cout << intervalSlots.String() << endl;
    //return 1;
}  


