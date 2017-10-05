
/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * 
 * https://github.com/Leofaber/AG_lm6
*/

#include <iostream>
#include <fstream>
#include <string.h>

 

#include <PilParams.h>
#include "BinEvaluator.h"
#include "LiMa.h"
#include "ExpRatioEvaluator.h"

using namespace std; 

const char* startString = {
"################################################################\n"
"###                   Task AG_lm6 v1.0.5 -               ###"
};

const char* endString = {
"### Task AG_lm6 exiting .................................... ###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "ctsT0", "Input T0 cts file name" },
    { PilString, "expT0", "Input T0 exp file name" },
    { PilString, "ctsT1", "Input T1 cts file name" },
    { PilString, "expT1", "Input T1 exp file name" },
    { PilString, "ctsT2", "Input T2 cts file name" },
    { PilString, "expT2", "Input T2 exp file name" },
    { PilReal, "l", "Longitude of GRB centroid (galactic)" },
    { PilReal, "b", "Latitude of GRB centroid (galactic)" },
    { PilReal, "radius", "Li&Ma radius of analysis" }, 
    { PilString, "normalize", "true if exp should be normalized, false otherwise" },
    { PilNone, "", "" }
};



int main(int argc, char *argv[])
{
    cout << startString << endl;
    
	float minTreshold = 0;
	float maxTreshold = 100;

    if(argc > 12) {
 		minTreshold = atof(argv[12]);
		maxTreshold = atof(argv[13]);
	}
		
	PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;
   
	const char * outfile = params["outfile"];
	const char *ctsT0FilePath = params["ctsT0"];
	const char *expT0FilePath = params["expT0"];
	const char *ctsT1FilePath = params["ctsT1"];
	const char *expT1FilePath = params["expT1"];
	const char *ctsT2FilePath = params["ctsT2"];
	const char *expT2FilePath = params["expT2"];
	double l = params["l"];
	double b = params["b"];
	double radius = params["radius"]; 
	const char *normalize = params["normalize"];

	bool doNormalization;
	
	if( strcmp(normalize, "true") == 0 )
		doNormalization = true;
	else
		doNormalization = false;
			
	
  /*	 FOR TESTING PURPOSES	
	const char *ctsT0FilePath = "../MAPFORTEST/T0.cts.gz";
	const char *expT0FilePath = "../MAPFORTEST/T0.exp.gz";
	const char *ctsT1FilePath = "../MAPFORTEST/T1.cts.gz";
	const char *expT1FilePath = "../MAPFORTEST/T1.exp.gz";
	const char *ctsT2FilePath = "../MAPFORTEST/T2.cts.gz";
	const char *expT2FilePath = "../MAPFORTEST/T2.exp.gz";
	double l = 45;
	double b = 30;
	double radius = 10;
	const char * outfile = "../outfile.txt";
	
*/

   	ofstream resText(outfile);
   	resText.setf(ios::fixed); 




	int statusCts = 0;
	int statusExp = 0;
	

    // EXPRATIOEVALUATOR OF EXPTO
	
	ExpRatioEvaluator expRatioT0(expT0FilePath); 
	double *expRatioArrayT0 = expRatioT0.computeExpRatioValues(l,b,doNormalization,minTreshold,maxTreshold); 
	if(expRatioArrayT0[0]!=-1) {
		cout << "ExpRatio evaluation of expT0: " << (int)round(expRatioArrayT0[0])<< endl;		
	}
	 
		
	// ANALYSIS OF SOURCE MAP T0
 	// Exp
	BinEvaluator expT0(expT0FilePath,l,b,radius);
	if(! expT0.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "expT0 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusExp = expT0.sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT0 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}
 	// Cts
	BinEvaluator ctsT0(ctsT0FilePath,l,b,radius);
	if(! ctsT0.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "ctsT0 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusCts = ctsT0.sumBin();
 	if(statusCts != 0)
	{
		fprintf(stderr,"ctsT0 Error: the radius exceeds the border of the .cts map\n");
		exit (EXIT_FAILURE);
	}

	if(statusCts == 0 && statusExp == 0) {
		resText << setprecision(1);
		resText << ctsT0.tmin << " " << ctsT0.tmax << " ";
		resText << setprecision(2);
		resText << ctsT0.binSum << " " << expT0.binSum << " ";
		resText << setprecision(10) << ctsT0.binSum / (double) expT0.binSum << " ";
		resText << setprecision(5); 
		resText << (int)round(expRatioArrayT0[0]) << " ";
	}
	
	
	
	
	
 
	// EXPRATIOEVALUATOR OF EXPT1

	ExpRatioEvaluator expRatioT1(expT1FilePath);
	double *expRatioArrayT1 = expRatioT1.computeExpRatioValues(l,b,doNormalization,minTreshold,maxTreshold); 
	if(expRatioArrayT1[0]!=-1) {
		cout << "ExpRatio evaluation of expT1: " << (int)round(expRatioArrayT1[0])<< endl;		
	}

	// ANALYSIS OF MAP T1

	BinEvaluator expT1(expT1FilePath,l,b,radius);
	if(! expT1.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "expT1 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusExp = expT1.sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT1 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}

 
	BinEvaluator ctsT1(ctsT1FilePath,l,b,radius);
	if(! ctsT1.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "ctsT1 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusCts = ctsT1.sumBin();
 	if(statusCts != 0)
	{
		fprintf(stderr,"ctsT1 Error: the radius exceeds the border of the .cts map\n");
		exit (EXIT_FAILURE);
	}

	if(statusCts == 0 && statusExp == 0) {
		resText << setprecision(1);
		resText << ctsT1.tmin << " " << ctsT1.tmax << " ";
		resText << setprecision(2);
		resText << ctsT1.binSum << " " << expT1.binSum << " ";
		resText << setprecision(5);
		resText << (int)round(expRatioArrayT1[0]) << " ";
	}
	
	
	// EXPRATIOEVALUATOR OF EXP T2

	ExpRatioEvaluator expRatioT2(expT2FilePath);
	double *expRatioArrayT2 = expRatioT2.computeExpRatioValues(l,b,doNormalization,minTreshold,maxTreshold); 
	if(expRatioArrayT2[0]!=-1) {
		cout << "ExpRatio evaluation of expT2: " << (int)round(expRatioArrayT2[0])<< endl;		
	}

	 
	// ANALYSIS OF MAP T2

	BinEvaluator expT2(expT2FilePath,l,b,radius);
	if(! expT2.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "expT2 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusExp = expT2.sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT2 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}

 
	BinEvaluator ctsT2(ctsT2FilePath,l,b,radius);
	if(! ctsT2.convertFitsDataToMatrix() )
	{
		fprintf( stderr, "ctsT2 convertFitsDataToMatrix() Error reading fits file\n");
		exit (EXIT_FAILURE);
	}

	statusCts = ctsT2.sumBin();
 	if(statusCts != 0)
	{
		fprintf(stderr,"ctsT2 Error: the radius exceeds the border of the .cts map\n");
		exit (EXIT_FAILURE);
	}

	if(statusCts == 0 && statusExp == 0) {
		resText << setprecision(1);
		resText << ctsT2.tmin << " " << ctsT2.tmax << " ";
		resText << setprecision(2);
		resText << ctsT2.binSum << " " << expT2.binSum << " ";
		resText << setprecision(5);
		resText << (int)round(expRatioArrayT2[0]) << " ";
		
	}
	
	


	//cout << ctsT0.binSum <<" "<<ctsT1.binSum <<" "<<ctsT2.binSum <<" "<<expT0.binSum <<" "<<expT1.binSum <<" "<<expT2.binSum<<endl;

	// LI&MA Analysis
	double S;
	cout << "\nLI&MA Analysis: " << endl;
	LiMa lm(ctsT0.binSum,ctsT1.binSum,ctsT2.binSum,expT0.binSum,expT1.binSum,expT2.binSum);
	
	if(expRatioArrayT0[0] != -1 && expRatioArrayT1[0] != -1 && expRatioArrayT2[0] != -1) { 
		S = lm.computeLiMiValue();
	}else{
		S=-1;
	}
	


	resText << lm.alpha << " " << std::setprecision(2)  << " off " << lm.bkg << " " << lm.expBgSum << " " << std::setprecision(10) << lm.bkg / (double) (lm.expBgSum)<< " " << S << endl;	//resText << SA << endl;
    
    resText.close();


    return 0;
}
