////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_checkMapValue
//       First release: 2015
//       Authors: Andrea Bulgarelli (INAF/OAS Bologna)
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

#include "PilParams.h"
#include "AgileMap.h"


using namespace std;

const PilDescription c_params[] = {
	{ PilString, "filename", "Map file name" },
	{ PilReal,   "l", "Longitude l in degrees (galactic)" },
	{ PilReal,   "b", "Latitude b in degrees (galactic)" },
	{ PilNone,   "",   "" }
	};


class CheckMapParams: public PilParams
{
public:
	CheckMapParams(): PilParams(c_params) {}
	void Print() { PilParams::Print(c_params); }
};


/**
static int DoPIL(
	int argc,
	char *argv[],
	char* filename,
	double& l,
	double& b)
{
int status = PILInit(argc, argv);
int numpar = 0;
status = status || PILGetNumParameters(&numpar);
GetPilParam("filename", filename, &status);
GetPilParam("l", &l, &status);
GetPilParam("b", &b, &status);
PILClose(status);
return status;
}
*/

int main(int argC, char* argV[])
{
CheckMapParams params;
if (!params.Load(argC, argV))
	return -1;

const char* filename = params["filename"];
double l = params["l"];
double b = params["b"];

AgileMap map;
if (!map.Read(filename)) {
	int row, col;
	bool inside = map.GetRowCol(l, b, &row, &col);
	// cout << "point (" << l << ", " << b << ") mapped to [" << row << ", " << col << "] " << (inside?"inside":"outside") << endl;
	if (inside) {
		cout << map(row, col) << endl;
		return 0;
		}
	else {
		cerr << "ERROR: point (" << l << ", " << b << ") falls outside the map" << endl;
		cout << -1 << endl;
		return -1;
		}
	}
else {
	cerr << "ERROR accessing file " << filename << endl;
	cout << -2 << endl;
	return -2;
	}
}
