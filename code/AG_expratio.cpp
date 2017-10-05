
/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * 
 * https://github.com/Leofaber/AG_expratio
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip> 
#include "ExpRatioEvaluator.h"
#include <cstring>


const char* startString = {
"################################################################\n"
"###                   Task AG_expratio v1.0.1 -              ###"
};

const char* endString = {
"### Task AG_expratio exiting ............................... ###\n"
"################################################################"
};

void printMatrix(double ** matrix, int rows, int cols){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j ++){
			cout << matrix[i][j] << " ";
		}
		cout << "\n";
	}
}


int main(int argc, char *argv[])
{ 	

	cout << startString << endl;
	
	if(argc < 6 || argc >8)
	{
		printf("\nAt least 5 arguments expected (+ 2 optional)\n   - The .exp file path\n   - The name of the output file\n   - Normalize boolean: true if exp-ratio must be computed on a normalized map, false otherwise\n   - The l coordinate\n    -The b coordinate\n\n(Optional)\n   - The minThreshold (default value = 0)\n   - The maxThreshold (default value = 100)\n");
		cout << endString << endl;		
		exit (EXIT_FAILURE);
	}

	else
	{	
        const char * imagePath = argv[1];
        const char * outfile = argv[2];
        const char *onNormalizedMap = argv[3];
		double l = atof(argv[4]);
		double b = atof(argv[5]);
		double minThreshold = 0;
		double maxThreshold = 100;
		if(argc == 8){
			minThreshold = (double) atof(argv[6]);
		    maxThreshold = (double) atof(argv[7]);
		}
		
		bool computeExpRatioOnNormalizedMap;
		if( strcmp(onNormalizedMap, "true") == 0 )
			computeExpRatioOnNormalizedMap = true;
		else
			computeExpRatioOnNormalizedMap = false;
		
		
		//cout << "imagePath: " << imagePath << endl;
		//cout << "outfile: " <<outfile << endl;
		//cout << "normalize value: "<<normalize << " doNormalization: " << doNormalization <<endl;
		//cout << "minThreshold: " << minThreshold << endl;
		//cout << "maxThreshold: " <<maxThreshold << endl;
		//cout << "l: " <<l << endl;
		//cout << "b: "<<b << endl;
		
		ofstream resText(outfile);
		resText.setf(ios::fixed); 

		ExpRatioEvaluator exp(imagePath);
		double *output = exp.computeExpRatioValues(l,b,computeExpRatioOnNormalizedMap,minThreshold,maxThreshold);

		resText << setprecision(5) << output[0];	// <<" "<< output[1] <<" "<< output[2] <<" "<< output[3]
		
		cout << "Created " <<outfile<< " log file."<<endl;
		resText.close();

		double ** expRatioMap = exp.createExpRatioPixelMap(computeExpRatioOnNormalizedMap,minThreshold,maxThreshold);
		
		//printMatrix(expRatioMap,exp.getRows(),exp.getCols());
	
		exp.writeMatrixDataInAgileMapFile(expRatioMap, exp.agileMap, "exp_ratio_norm.exp");	
		
		cout << endString << endl;
		return 0;

	}
}
		
		
