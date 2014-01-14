


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
