
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



#include <iomanip> 
#include <cstring>

#include "ExpRatioEvaluator.h"



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



	// DEFINIZIONE PARAMETRI DI DEFAULT
	double minThreshold = 120;
	double maxThreshold = 140;
	int squareSize = 10;

	// CONTROLLO NUMERO PARAMETRI (TOO FEW, TOO MUCH)	
	if(argc < 6 || argc > 9)
	{
		printf("\nAt least 5 arguments expected (+ 3 optional)\n   - The .exp file path\n   - The name of the output file\n   - Normalize boolean: true if exp-ratio must be computed on a normalized map, false otherwise\n   - The l coordinate\n    -The b coordinate\n\n(Optional)\n   - The minThreshold (default value = 120)\n   - The maxThreshold (default value = 140)\n   - The square size (default value = 10)\n");
		cout << endString << endl;		
		exit (EXIT_FAILURE);
	}

	else
	{	

		// PARAMETRI OBBLIGATORI
        const char * imagePath = argv[1];
        const char * outfile = argv[2];
        const char *onNormalizedMap = argv[3];
		double l = atof(argv[4]);
		double b = atof(argv[5]);

		// CONTROLLO PARAMETRI OPZIONALI
		if(argc == 7)
		{	
			if(((string)argv[6])!="d")
				minThreshold = atof(argv[6]);
		}
		else if(argc == 8)
		{
			if(((string)argv[6])!="d")
				minThreshold = atof(argv[6]);
			if(((string)argv[7])!="d")
				maxThreshold = atof(argv[7]);
		}
		else if(argc == 9)
		{
			if(((string)argv[6])!="d")
				minThreshold = atof(argv[6]);
			if(((string)argv[7])!="d")
				maxThreshold = atof(argv[7]);
			if(((string)argv[8])!="d")
				squareSize = atof(argv[8]);
		}

				
		
		
		cout << "MinThreshold: " << minThreshold << endl;
		cout << "MaxThreshold: " << maxThreshold << endl;
		cout << "squareSize: " << squareSize << endl;
		bool computeExpRatioOnNormalizedMap;
		if( strcmp(onNormalizedMap, "true") == 0 ){
			computeExpRatioOnNormalizedMap = true;
			
		}
		else{
			computeExpRatioOnNormalizedMap = false;
			
		}
		
		
		ofstream resText(outfile);
		resText.setf(ios::fixed); 


		ExpRatioEvaluator exp(imagePath,computeExpRatioOnNormalizedMap,minThreshold,maxThreshold,squareSize);
		
					

		double expRatio = exp.computeExpRatioValues(l,b);

		resText << setprecision(5) << expRatio; 
		
		cout << "Created " <<outfile<< " log file."<<endl;
		resText.close();
		
		cout << endString << endl;
		return 0;

	}
}
		
		
