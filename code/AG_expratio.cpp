
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
"###                   Task AG_expratio v1.0.4 -              ###"
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


	// DEFINIZIONE PARAMETRI DI DEFAULT ---------------------------------------------


	double minThreshold = 120;
	double maxThreshold = 140;
	int squareSize = 20;
	const char * isAlreadyNormalized = "false";
	const char * createExpNormalizedMap = "false";
	const char * createExpRatioMap = "false";





	// CONTROLLO NUMERO PARAMETRI (TOO FEW, TOO MUCH) -------------------------------

	
	if(argc < 6 || argc > 11)
	{
		printf("\nAt least 5 arguments expected (+ 5 optional)\n   - The name of the output file\n   - The .exp file path\n   - The isExpMapNormalized boolean: if 'true' we assert that the exp map given in input is already normalized.\n   - The l coordinate\n   - The b coordinate\n\n(Optional)\n   - The createExpNormalizedMap: if 'true' the ExpNormalizedMap will be writed on file. (default value = false)\n   - The createExpRatioMap boolean: if 'true' the ExpRatioMap will be writed on file. (default value = false)\n   - The minThreshold (default value = 120)\n   - The maxThreshold (default value = 140)\n   - The square size (default value = 20).\n");
		cout << endString << endl;		
		exit (EXIT_FAILURE);
	}

	else
	{	

		// PARAMETRI OBBLIGATORI ---------------------------------------------------
        const char * outfile = argv[1];
        const char * imagePath = argv[2];
        isAlreadyNormalized = argv[3];
		double l = atof(argv[4]);
		double b = atof(argv[5]);

        
		

		// PARAMETRI OPZIONALI -------------------------------------------
		if(argc == 7)
		{	
			if(((string)argv[6])!="d")
				createExpNormalizedMap = argv[6];
		}
		else if(argc == 8)
		{
			if(((string)argv[6])!="d")
				createExpNormalizedMap = argv[6];
			if(((string)argv[7])!="d")
				createExpRatioMap = argv[7];
		}
		else if(argc == 9)
		{	
			if(((string)argv[6])!="d")
				createExpNormalizedMap = argv[6];
			if(((string)argv[7])!="d")
				createExpRatioMap = argv[7];
			if(((string)argv[8])!="d")
				minThreshold = atof(argv[8]);
			}
		
		else if(argc == 10)
		{
			if(((string)argv[6])!="d")
				createExpNormalizedMap = argv[6];
			if(((string)argv[7])!="d")
				createExpRatioMap = argv[7];
			if(((string)argv[8])!="d")
				minThreshold = atof(argv[8]);
			if(((string)argv[9])!="d")
				maxThreshold = atof(argv[9]);
	
		}
		else if(argc == 11)
		{
			if(((string)argv[6])!="d")
				createExpNormalizedMap = argv[6];
			if(((string)argv[7])!="d")
				createExpRatioMap = argv[7];
			if(((string)argv[8])!="d")
				minThreshold = atof(argv[8]);
			if(((string)argv[9])!="d")
				maxThreshold = atof(argv[9]);
			if(((string)argv[10])!="d")
				squareSize = atof(argv[10]);
		}
		
		bool isAlreadyNormalizedBool;
		if( strcmp(isAlreadyNormalized, "true") == 0 ){
			isAlreadyNormalizedBool = true;
		}
		else{
			isAlreadyNormalizedBool = false;
		}				
		
		bool createExpNormalizedMapBool;
		if( strcmp(createExpNormalizedMap, "true") == 0 ){
			createExpNormalizedMapBool = true;
		}
		else{
			createExpNormalizedMapBool = false;
		}
		
		bool createExpRatioMapBool;
		if( strcmp(createExpRatioMap, "true") == 0 ){
			createExpRatioMapBool = true;
		}
		else{
			createExpRatioMapBool = false;
		}
		
		



		
		// PRINT INPUT PARAMETERS -------------------------------------


		cout << "\noutfile: " << outfile << endl;
		cout << "imagePath: " << imagePath << endl;
		cout << "isAlreadyNormalized: " << isAlreadyNormalized << endl;
		cout << "l: " << l << endl;
		cout << "b: " << b << endl;
		cout << "createExpNormalizedMap: " << createExpNormalizedMap << endl;
		cout << "createExpRatioMap: " << createExpRatioMap << endl;
		cout << "MinThreshold: " << minThreshold << endl;
		cout << "MaxThreshold: " << maxThreshold << endl;
		cout << "squareSize: " << squareSize << "\n" <<endl;
		



		// CORE LOGIC -----------------------------------------------


		ExpRatioEvaluator exp(imagePath, isAlreadyNormalizedBool, createExpNormalizedMapBool, createExpRatioMapBool, minThreshold, maxThreshold, squareSize);
					
		double expRatio = exp.computeExpRatioValues(l,b);







		// OUTPUT ----------------------------------------------------
		
		
		ofstream resText(outfile);
		resText.setf(ios::fixed); 

		resText << setprecision(5) << expRatio; 
		
		cout << "Created " <<outfile<< " log file."<<endl;
		resText.close();
		
		cout << endString << endl;
		return 0;

	}
}
		
		
