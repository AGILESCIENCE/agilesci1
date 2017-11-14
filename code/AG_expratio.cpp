
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


#include <fstream>
#include <iomanip> 
#include <cstring>
#include <PilParams.h>

#include "ExpRatioEvaluator.h"


const char* startString = {
"################################################################\n"
"###                   Task AG_expratio v1.0.5 -              ###"
};

const char* endString = {
"### Task AG_expratio exiting ............................... ###\n"
"################################################################"
};

const PilDescription paramsDescr[] = {
	{ PilString, "outfile","Output file name"},
	{ PilString, "imagePath", "Input expMap file name"},
	{ PilBool, "isExpMapNormalized","If 'yes' (or 'y') you assert that the exp-map is already normalized. Insert 'None' (or 'no') instead"},
	{ PilReal, "l", "Longitude of map center (galactic)" },
    	{ PilReal, "b", "Latitude of map center (galactic)" },        
    	{ PilBool, "createExpNormalizedMap", "If 'yes' (or 'y') the normalized exp map will be written on file"},
    	{ PilBool, "createExpRatioMap", "If true the exp-ratio map will be written on file"},
    	{ PilReal, "minThreshold", "The lower bound for the threshold level in exp-ratio evaluation"},
    	{ PilReal, "maxThreshold", "The upper bound for the threshold level in exp-ratio evaluation"},
    	{ PilReal, "squareSize", "The degree dimension of the exp ratio evaluation area's edge"},
   	{ PilNone, "", "" }
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
	
	
	PilParams params(paramsDescr);
   	if (!params.Load(argc, argv))
        	return EXIT_FAILURE;
        
        
    	// PARAMETRI OBBLIGATORI ---------------------------------------------------
    	const char *outfile		= params["outfile"];
    	const char *imagePath 		= params["imagePath"];	
    	bool isExpMapNormalized 	= params["isExpMapNormalized"];	
	double l 			= params["l"];	
	double b 			= params["b"];	
	bool createExpNormalizedMap 	= params["createExpNormalizedMap"];
	bool createExpRatioMap 		= params["createExpRatioMap"];
	double minThreshold 		= params["minThreshold"]; 
	double maxThreshold 		= params["maxThreshold"]; 
	double squareSize 		= params["squareSize"];  
	 

	// PRINT INPUT PARAMETERS -------------------------------------


	cout << "\noutfile: " 			<< outfile << endl;
	cout << "imagePath: " 			<< imagePath << endl;
	cout << "isExpMapNormalized: " 		<< isExpMapNormalized << endl;
	cout << "l: "				<< l << endl;
	cout << "b: " 				<< b << endl;
	cout << "createExpNormalizedMap: " 	<< createExpNormalizedMap << endl;
	cout << "createExpRatioMap: " 		<< createExpRatioMap << endl;
	cout << "MinThreshold: " 		<< minThreshold << endl;
	cout << "MaxThreshold: "		<< maxThreshold << endl;
	cout << "squareSize: " 			<< squareSize << "\n" <<endl;
	



	// CORE LOGIC -----------------------------------------------


	ExpRatioEvaluator exp(	imagePath,
				isExpMapNormalized,
				createExpNormalizedMap, 
				createExpRatioMap,
				minThreshold, 
				maxThreshold,
				squareSize);
				
	double expRatio = exp.computeExpRatioValues(l,b);


	// OUTPUT ----------------------------------------------------
	
 
	ofstream resText(outfile);
	resText.setf(ios::fixed); 

	resText << setprecision(5) << expRatio; 
	
	cout << "Created " << outfile << " log file."<<endl;
	resText.close();
	
	cout << endString << endl;
	return 0;

	
}
		
		
