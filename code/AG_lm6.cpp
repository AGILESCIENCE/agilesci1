
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
"###                   Task AG_lm6 v1.0.8 -               ###"
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
    { PilString, "isExpMapsNormalized", "insert true if T0,T1,T2 exp maps are already normalized, insert false otherwise" },
    { PilReal, "l", "Longitude of GRB centroid (galactic)" },
    { PilReal, "b", "Latitude of GRB centroid (galactic)" },
    { PilReal, "radius", "Li&Ma radius of analysis" },
    { PilNone, "", "" }
};



int main(int argc, char *argv[])
{
    cout << startString << endl;

	// PARAMETRI OPZIONALI - VALORI DI DEFAULT
	const char * binSumOnNormalizedMap = "true";
    	bool binSumOnNormalizedMapBool = true;

	const char * createExpNormalizedMap = "false";
	bool createExpNormalizedMapBool = false;

	const char * createExpRatioMap = "false";
	bool createExpRatioMapBool = false;

	float minThreshold = 120;
	float maxThreshold = 140;
	int squareSize = 20;
	
	
	 
	

	// CONTROLLO NUMERO PARAMETRI (TOO FEW, TOO MUCH)	
	if(argc < 12 || argc > 18)
	{
		printf("\nAt least 11 arguments expected (+ 6 optional)\n   - The name of the output file\n   - The Input T0 cts file name\n   - The Input T0 exp file name\n   - The Input T1 cts file name\n   - The Input T1 exp file name\n   - The Input T2 cts file name\n   - The Input T2 exp file name\n   - insert true if T0,T1,T2 exp maps are already normalized, insert false otherwise\n   - The l coordinate\n   - The b coordinate\n   - The radius size\n\n(Optional)\n   - compute bin sum on normalized maps (default =true)\n   - createExpNormalizedMap: insert true to write the normalized maps for T0,T1 and T2 on file (default = false)\n   - createExpRatioMap: insert true to write the exp-ratio maps for T0,T1 and T2 on file (default = false)\n   - The minThreshold (default value = 120)\n   - The maxThreshold (default value = 140)\n   - The square size (default value = 20)\n\n ");
		cout << endString << endl;		
		exit (EXIT_FAILURE);
	}

	PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;
	

	// PARAMETRI OBBLIGATORI
	const char *outfile = params["outfile"];

	const char *ctsT0FilePath = params["ctsT0"];
	const char *expT0FilePath = params["expT0"];

	const char *ctsT1FilePath = params["ctsT1"];
	const char *expT1FilePath = params["expT1"];

	const char *ctsT2FilePath = params["ctsT2"];
	const char *expT2FilePath = params["expT2"];
	
	const char * isExpMapsNormalizedChar = params["isExpMapsNormalized"];
	bool isExpMapsNormalizedBool = false;
	if( strcmp(isExpMapsNormalizedChar, "true") == 0 )
		isExpMapsNormalizedBool = true;

	double l = params["l"];
	double b = params["b"];
	double radius = params["radius"]; 

	
	
	
	


	
	
		

	// PARAMETRI OPZIONALI
    if(argc == 13) 	
    {
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		
    }
	if(argc == 14) 
	{
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		if(((string)argv[13])!="d")
			createExpNormalizedMap = argv[13];
	}
	if(argc == 15) 
	{
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		if(((string)argv[13])!="d")
			createExpNormalizedMap = argv[13];
		if(((string)argv[14])!="d")		
			createExpRatioMap = argv[14];
		
	}
	if(argc == 16) 
	{
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		if(((string)argv[13])!="d")
			createExpNormalizedMap = argv[13];
		if(((string)argv[14])!="d")		
			createExpRatioMap = argv[14];
		if(((string)argv[15])!="d")
			minThreshold = atof(argv[15]);
	}	
	if(argc == 17) 
	{
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		if(((string)argv[13])!="d")
			createExpNormalizedMap = argv[13];
		if(((string)argv[14])!="d")		
			createExpRatioMap = argv[14];
		if(((string)argv[15])!="d")
			minThreshold = atof(argv[15]);
		if(((string)argv[16])!="d")		
			maxThreshold = atof(argv[16]);
		
	}	
	if(argc == 18) 
	{
		if(((string)argv[12])!="d")
			binSumOnNormalizedMap = argv[12];
		if(((string)argv[13])!="d")
			createExpNormalizedMap = argv[13];
		if(((string)argv[14])!="d")		
			createExpRatioMap = argv[14];
		if(((string)argv[15])!="d")
			minThreshold = atof(argv[15]);
		if(((string)argv[16])!="d")		
			maxThreshold = atof(argv[16]);
		if(((string)argv[17])!="d")		
			squareSize = atoi(argv[17]);
		
	}
		
    

	if( strcmp(binSumOnNormalizedMap, "true") == 0 )
		binSumOnNormalizedMapBool = true;
	else
		binSumOnNormalizedMapBool = false;

	if( strcmp(createExpNormalizedMap, "true") == 0 )
		createExpNormalizedMapBool = true;

	if( strcmp(createExpRatioMap, "true") == 0 )
		createExpRatioMapBool = true;



	 
	// INPUT PARAMETERS ---------------------------------------------------

	cout << "\noutfile: " << outfile << endl;
	cout << "ctsT0FilePath: " << ctsT0FilePath << endl;
	cout << "expT0FilePath: " << expT0FilePath << endl;
	cout << "ctsT1FilePath: " << ctsT1FilePath << endl;
	cout << "expT1FilePath: " << expT1FilePath << endl;
	cout << "ctsT2FilePath: " << ctsT2FilePath << endl;
	cout << "expT2FilePath: " << expT2FilePath << endl;
	cout << "isExpMapsNormalizedBool: " << isExpMapsNormalizedBool << endl;
	cout << "l: " << l << endl;
	cout << "b: " << b << endl;
	cout << "radius: " << radius << endl;
	cout << "binSumOnNormalizedMapBool: " << binSumOnNormalizedMapBool << endl;
	cout << "createExpNormalizedMapBool: " << createExpNormalizedMapBool << endl;
	cout << "createExpRatioMapBool: " << createExpRatioMapBool << endl;
	cout << "minThreshold: " << minThreshold << endl;
	cout << "maxThreshold: " << maxThreshold << endl;
	cout << "squareSize: " << squareSize << "\n"<<endl;
	

 
	ofstream resText(outfile);
   	resText.setf(ios::fixed); 




	int statusCts = 0;
	int statusExp = 0;
	
	

    // EXPRATIOEVALUATOR OF EXPT0
	
	ExpRatioEvaluator expRatioT0(expT0FilePath, isExpMapsNormalizedBool, createExpNormalizedMapBool, createExpRatioMapBool, minThreshold, maxThreshold, squareSize);
	double expRatioValueT0 = expRatioT0.computeExpRatioValues(l,b); 
	if(expRatioValueT0!=-1) { 
		cout << "ExpRatio evaluation of expT0: " << (int)round(expRatioValueT0)<< endl;		
	}
	 
		
	// ANALYSIS OF SOURCE MAP T0
 	// Exp
 	BinEvaluator * beT0;
	if(binSumOnNormalizedMapBool)
		beT0 = new BinEvaluator(expT0FilePath,expRatioT0.getNormalizedMap(),l,b,radius);
	else
		beT0 = new BinEvaluator(expT0FilePath,expRatioT0.getImage(),l,b,radius);
	

	
	
	statusExp = beT0->sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT1 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}
 
	
 	// Cts
	BinEvaluator ctsT0(ctsT0FilePath,l,b,radius);
	
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
		resText << (int) ctsT0.binSum << " " << beT0->binSum << " ";
		resText << setprecision(10) << ctsT0.binSum / (double) beT0->binSum << " ";
		resText << setprecision(5); 
		resText << (int)round(expRatioValueT0) << " ";		
	}
	
	
	
	
	
 
	// EXPRATIOEVALUATOR OF EXPT1

	ExpRatioEvaluator expRatioT1(expT1FilePath, isExpMapsNormalizedBool, createExpNormalizedMapBool, createExpRatioMapBool, minThreshold, maxThreshold, squareSize);
	double expRatioValueT1 = expRatioT1.computeExpRatioValues(l,b); 
	if(expRatioValueT1!=-1) {	
		cout << "ExpRatio evaluation of expT1: " << (int)round(expRatioValueT1)<< endl;				
	}

	// ANALYSIS OF MAP T1
	BinEvaluator * beT1;
	if(binSumOnNormalizedMapBool) 
		beT1 = new BinEvaluator(expT1FilePath,expRatioT1.getNormalizedMap(),l,b,radius);
	else
		beT1 = new BinEvaluator(expT1FilePath,expRatioT1.getImage(),l,b,radius);
	

	statusExp = beT1->sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT1 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}
	

	

 
	BinEvaluator ctsT1(ctsT1FilePath,l,b,radius);
	

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
		resText << (int) ctsT1.binSum << " " << beT1->binSum << " ";
		resText << setprecision(5);
		resText << (int)round(expRatioValueT1) << " ";	
	}
	
	
	// EXPRATIOEVALUATOR OF EXP T2

	ExpRatioEvaluator expRatioT2(expT2FilePath, isExpMapsNormalizedBool, createExpNormalizedMapBool, createExpRatioMapBool, minThreshold, maxThreshold, squareSize);
	double expRatioValueT2 = expRatioT2.computeExpRatioValues(l,b); 
	if(expRatioValueT2!=-1) {
		cout << "ExpRatio evaluation of expT2: " << (int)round(expRatioValueT2)<< endl;			
	}

	 
	// ANALYSIS OF MAP T2
	BinEvaluator * beT2;
	if(binSumOnNormalizedMapBool) 
		beT2 = new BinEvaluator(expT2FilePath,expRatioT2.getNormalizedMap(),l,b,radius);
	else
		beT2 = new BinEvaluator(expT2FilePath,expRatioT2.getImage(),l,b,radius);
	

	statusExp = beT2->sumBin();
	if(statusExp != 0)
	{
		fprintf(stderr,"expT2 Error: the radius exceeds the border of the .exp map\n");
		exit (EXIT_FAILURE);
	}
	

 
	BinEvaluator ctsT2(ctsT2FilePath,l,b,radius);
	

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
		resText << (int) ctsT2.binSum << " " << beT2->binSum << " ";
		resText << setprecision(5);
		resText << (int)round(expRatioValueT2) << " ";		
		
	}
	
	


	//cout << ctsT0.binSum <<" "<<ctsT1.binSum <<" "<<ctsT2.binSum <<" "<<expT0.binSum <<" "<<expT1.binSum <<" "<<expT2.binSum<<endl;

	// LI&MA Analysis
	double S;
	cout << "\nLI&MA Analysis: " << endl;
	LiMa lm(ctsT0.binSum,ctsT1.binSum,ctsT2.binSum,beT0->binSum,beT1->binSum,beT2->binSum);
	
	//if(expRatioArrayT0 != -1 && expRatioArrayT1 != -1 && expRatioArrayT2 != -1) { 		//elimintao [0]
		S = lm.computeLiMiValue();
	/*}else{
		S=-1;
	}*/
	


	resText << lm.alpha << " " << std::setprecision(2)  << " off " << lm.bkg << " " << lm.expBgSum << " " << std::setprecision(10) << lm.bkg / (double) (lm.expBgSum)<< " " << S << endl;	//resText << SA << endl;
    
    resText.close();

	cout << endString << endl;

    return 0;
  
}
