////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_lm6
//       Author: Andrea Bulgarelli, Leonardo Baroncelli, Giancarlo Zollino
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
#include <fstream>
#include <string.h>
#include <PilParams.h>

#include "BinEvaluator.h"
#include "LiMa.h"
#include "ExpRatioEvaluator.h"

using namespace std;

const char* startString = {
"################################################################\n"
"###              AG_lm6 B25 v1.0.9 - L.B. G.Z.               ###\n"
"################################################################\n"
};

const char* endString = {
"################################################################\n"
"###  AG_lm6 exiting B25 .................................... ###\n"
"################################################################\n"
};

const PilDescription paramsDescr[] = {
    { PilString, "outfile", "Output file name" },
    { PilString, "ctsT0", "Input T0 cts file name" },
    { PilString, "expT0", "Input T0 exp file name" },
    { PilString, "ctsT1", "Input T1 cts file name" },
    { PilString, "expT1", "Input T1 exp file name" },
    { PilString, "ctsT2", "Input T2 cts file name" },
    { PilString, "expT2", "Input T2 exp file name" },
    { PilBool, "isExpMapsNormalized", "insert true if T0,T1,T2 exp maps are already normalized, insert false otherwise" },
    { PilReal, "l", "Longitude of GRB centroid (galactic)" },
    { PilReal, "b", "Latitude of GRB centroid (galactic)" },
    { PilReal, "radius", "Li&Ma radius of analysis" },
    { PilBool, "binSumOnNormalizedMap","compute bin sum on normalized maps (default =true)"},
    { PilBool, "createExpNormalizedMap","If 'yes' (or 'y') the normalized exp maps will be written on file"},
    { PilBool, "createExpRatioMap", "If 'yes' (or 'y') the exp-ratio maps will be written on file"},
    { PilReal, "minThreshold", "The lower bound for the threshold level in exp-ratio evaluation"},
    { PilReal, "maxThreshold", "The upper bound for the threshold level in exp-ratio evaluation"},
    { PilReal, "squareSize", "The degree dimension of the exp ratio evaluation area's edge"},
    { PilNone, "", "" }
};



int main(int argc, char *argv[])
{
  cout << startString << endl;

	PilParams params(paramsDescr);
    if (!params.Load(argc, argv))
        return EXIT_FAILURE;


	// PARAMETRI OBBLIGATORI ---------------------------------------------------

	const char *outfile = params["outfile"];

	const char *ctsT0FilePath = params["ctsT0"];
	const char *expT0FilePath = params["expT0"];

	const char *ctsT1FilePath = params["ctsT1"];
	const char *expT1FilePath = params["expT1"];

	const char *ctsT2FilePath = params["ctsT2"];
	const char *expT2FilePath = params["expT2"];

	bool isExpMapsNormalized = params["isExpMapsNormalized"];

	double l = params["l"];
	double b = params["b"];
	double radius = params["radius"];

	bool binSumOnNormalizedMap = params["binSumOnNormalizedMap"];
 	bool createExpNormalizedMap = params["createExpNormalizedMap"];
	bool createExpRatioMap = params["createExpRatioMap"];
	double minThreshold = params["minThreshold"];
	double maxThreshold = params["maxThreshold"];
	double squareSize = params["squareSize"];

  params.Print();

	ofstream resText(outfile);
   	resText.setf(ios::fixed);

	int statusCts = 0;
	int statusExp = 0;


  /////////////////////////////////////////////////////////////////////////////
  //
  //  Exp-ratio evaluation of expT0

  	ExpRatioEvaluator expRatioT0(expT0FilePath, isExpMapsNormalized, createExpNormalizedMap, createExpRatioMap, minThreshold, maxThreshold, squareSize);
  	double expRatioValueT0 = expRatioT0.computeExpRatioValues(l,b);
  	if(expRatioValueT0!=-1)
  		cout << "ExpRatio evaluation of expT0: " << (int)round(expRatioValueT0) << endl;



  /////////////////////////////////////////////////////////////////////////////
  //
  //  Bin Evaluation on ctsT0 and expT0

      // Exp
   	BinEvaluator * beT0;

  	binSumOnNormalizedMap ? beT0 = new BinEvaluator(expT0FilePath,expRatioT0.getNormalizedMap(),l,b,radius) : beT0 = new BinEvaluator(expT0FilePath,expRatioT0.getImage(),l,b,radius);

  	beT0->sumBin();

   	// Cts
  	BinEvaluator ctsT0(ctsT0FilePath,l,b,radius);

    ctsT0.sumBin();





  	resText << setprecision(1);
  	resText << ctsT0.tmin << " " << ctsT0.tmax << " ";
  	resText << setprecision(2);
  	resText << (int) ctsT0.binSum << " " << beT0->binSum << " ";
  	resText << setprecision(10) << ctsT0.binSum / (double) beT0->binSum << " ";
  	resText << setprecision(5);
  	resText << (int)round(expRatioValueT0) << " ";


  /////////////////////////////////////////////////////////////////////////////
  //
  //  Exp-ratio evaluation of expT1

  	ExpRatioEvaluator expRatioT1(expT1FilePath, isExpMapsNormalized, createExpNormalizedMap, createExpRatioMap, minThreshold, maxThreshold, squareSize);
  	double expRatioValueT1 = expRatioT1.computeExpRatioValues(l,b);
  	if(expRatioValueT1!=-1)
  		cout << "ExpRatio evaluation of expT1: " << (int)round(expRatioValueT1)<< endl;


  /////////////////////////////////////////////////////////////////////////////
  //
  //  Bin Evaluation on ctsT0 and expT0

      // Exp
  	BinEvaluator * beT1;
  	binSumOnNormalizedMap ? beT1 = new BinEvaluator(expT1FilePath,expRatioT1.getNormalizedMap(),l,b,radius) : beT1 = new BinEvaluator(expT1FilePath,expRatioT1.getImage(),l,b,radius);
      beT1->sumBin();

      // Cts
  	BinEvaluator ctsT1(ctsT1FilePath,l,b,radius);
  	ctsT1.sumBin();

  	resText << setprecision(1);
  	resText << ctsT1.tmin << " " << ctsT1.tmax << " ";
  	resText << setprecision(2);
  	resText << (int) ctsT1.binSum << " " << beT1->binSum << " ";
  	resText << setprecision(5);
  	resText << (int)round(expRatioValueT1) << " ";




  /////////////////////////////////////////////////////////////////////////////
  //
  //  Exp-ratio evaluation of expT2

  	ExpRatioEvaluator expRatioT2(expT2FilePath, isExpMapsNormalized, createExpNormalizedMap, createExpRatioMap, minThreshold, maxThreshold, squareSize);
  	double expRatioValueT2 = expRatioT2.computeExpRatioValues(l,b);
  	if(expRatioValueT2!=-1)
  		cout << "ExpRatio evaluation of expT2: " << (int)round(expRatioValueT2)<< endl;

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Bin Evaluation on ctsT0 and expT0

      // Exp
  	BinEvaluator * beT2;
  	binSumOnNormalizedMap ? beT2 = new BinEvaluator(expT2FilePath,expRatioT2.getNormalizedMap(),l,b,radius) : beT2 = new BinEvaluator(expT2FilePath,expRatioT2.getImage(),l,b,radius);
      beT2->sumBin();

      // Cts
  	BinEvaluator ctsT2(ctsT2FilePath,l,b,radius);
      ctsT2.sumBin();

    resText << setprecision(1);
  	resText << ctsT2.tmin << " " << ctsT2.tmax << " ";
  	resText << setprecision(2);
  	resText << (int) ctsT2.binSum << " " << beT2->binSum << " ";
  	resText << setprecision(5);
  	resText << (int)round(expRatioValueT2) << " ";


  /////////////////////////////////////////////////////////////////////////////
  //
  //  Li & Ma analysis
    cout << "\nLI&MA Analysis: " << endl;

  	LiMa lm(ctsT0.binSum,ctsT1.binSum,ctsT2.binSum,beT0->binSum,beT1->binSum,beT2->binSum);

    double S = lm.computeLiMiValue();



  /////////////////////////////////////////////////////////////////////////////
  //
  //  Output

  	resText << lm.alpha << " " << std::setprecision(2)  << " off " << lm.bkg << " " << lm.expBgSum << " " << std::setprecision(10) << lm.bkg / (double) (lm.expBgSum)<< " " << S << endl;	//resText << SA << endl;

    resText.close();

  	cout << endString << endl;


    return 0;

}
