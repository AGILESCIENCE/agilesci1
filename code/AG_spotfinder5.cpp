





#define NEW_VERSION

#ifdef NEW_VERSION

#include <iostream>
#include <fstream>


#include <cstdlib>
#include "AgileMap.h"
#include "PlotCts2D3.h"

using namespace std;

#define AlikeMap AgileMap
#define GetNrows Rows
#define GetNcols Cols
#define AlikeSphdistDeg SphDistDeg

#else

#include <cstdlib>
#include "AlikeMap4.h"
#include "PlotCts2D3.h"

#endif



int main(int argc, char* argv[])
{

if (argc == 4) { // si stampa in stdout le coordinate dei pixels
/*	cout << "0) input file fits" << endl;
	cout << "1) input file txt " << endl;
	cout << "2) output file " << endl;*/
	const char* file1 = argv[1];
	const char* file2 = argv[2];
	const char* file3 = argv[3];
	AlikeMap map(file1);
	double lcenter = map.GetMapCenterL();
	double bcenter = map.GetMapCenterB();
	ifstream inFile2(file2);
	if (!inFile2) {
		cout << "Unable to open file (2)";
		return -1;	// terminate with error
		}

	ofstream asciiFile(file3);

	double L, B, TS, COUNTS, GAL, ISO;
	while (inFile2 >> L >> B >> TS >> COUNTS >> GAL >> ISO) {
		int px, py;	
		bool insideMap = map.GetRowCol(L, B, &px, &py);	/// What to do if outside the map?
		double dist = AlikeSphdistDeg(L, B, lcenter, bcenter);
		asciiFile << L << " " << B << " " << TS << " " << COUNTS << " "
					<< GAL << " " << ISO << " " << px << " " << py << " " << dist << endl;	
		}
	return 0;
	}

if (argc == 6) {
	cout << "0) input file fits" << endl;
	cout << "1) output file txt " << endl;
	cout << "2) radious of extraction " << endl;
	cout << "3) step of bins " << endl;
	cout << "4) type of analysis (1 or 2) " << endl;
	const char* file1 = argv[1];
	const char* file2 = argv[2];
	float radious = atof(argv[3]);
	float stepbin = atof(argv[4]);	/// zzz Questo sembra dover essere un intero
	int typeanal = atoi(argv[5]);

	AlikeMap map(file1);

	double lcenter = map.GetMapCenterL();
	double bcenter = map.GetMapCenterB();

	std::ofstream asciiFile(file2);
	int index = 10000;

	int M = map.GetNrows();
	int N = map.GetNcols();
	for (int i=0; i<M; i+=stepbin)
		for (int j=0; j<N; j+=stepbin) {
			double l = map.l(i, j);
			double b = map.b(i, j);
			double dist = AlikeSphdistDeg(l, b, lcenter, bcenter);
			if (dist<=radious) 
				asciiFile << "0.00000e-08 " << l << " " << b << " " << " 2.1 " << typeanal << " 2.0 " << index++ << endl;
			}
	return 0;
	}

	
if (argc == 1) {	/// zzz ???
	cout << "0) input file" << endl;
	cout << "1) binsize of the input " << endl;
	cout << "2) smoothing (n bins) " << endl;
	cout << "3) max number of connected regions to found " << endl;
	cout << "4) output files " << endl;
	cout << "5) algorithm type (0 original, 1 new with baricenter calculation, 2 new without baricenter calculation) " << endl;
	cout << "6) remove spot too near (radious, if 0, dont remove)" << endl;
	cout << "7) (optional) sky segmentation (0 - all sky, 1 - b >= 10, 2 - |b|<10, 3 - b <= -10)" << endl;
	cout << "8) (optional) shift the result coordinate to north - true or false (if true, b = b + bin size)" << endl;
	cout << "9) (optional) remove sources outside radius r from the center of the map (default 0, don't remove)" << endl;
	cout << "10) (optional) exposure file name (default "", don't use)" << endl;
	cout << "11) (optional) min exposure (default 200)" << endl;
	return -1;
	}
cout << "Number of argument " << argc << endl;

/// Mandatory arguments
const char* file = argv[1];
float binsize = atof(argv[2]);
int smoothing = atoi(argv[3]); 
int NCONREGSEARCH = atoi(argv[4]);
const char* out = argv[5];
int type = atoi(argv[6]);
double radiousremove = atof(argv[7]);

/// Optional arguments and their default
int skysegmentation = argc>8 ? atoi(argv[8]) : 0;
if (skysegmentation<0 || skysegmentation>3) {
	cerr << "Option 8 sky segmentation can be 0, 1, 2 or 3" << endl;
	return -1;
	}
int shiftnorth_b = argc>9 ? atoi(argv[9]) : 0;
double sourceradiusremove = argc>10 ? atof(argv[10]) : 0;
const char* expfile = argc>11 ? argv[11] : "";
double minexp = argc>12 ? atof(argv[12]) : 200;

	
// 	TApplication theApp("A", &argc, argv);
// 	theApp.SetReturnFromRun(kTRUE); 

if (type == 0)
	PlotCts2D(shiftnorth_b, skysegmentation, file, binsize, sourceradiusremove, smoothing, NCONREGSEARCH, out);
if (type == 1 || type == 2)
	PlotCts2D_FINE(shiftnorth_b, skysegmentation, file, binsize, sourceradiusremove, smoothing, NCONREGSEARCH, out, radiousremove, expfile, minexp);
// 	theApp.Run();  

}

