


#include <iostream>
#include <cstdlib>

#include "AgileMap.h"
#include "HealpixMap.h"


using std::cerr;
using std::cout;
using std::endl;




static int Help()
{
cout << "Usage: cts2healpix [+/-] <input_name> <output_name> [<nside>=512] [<rowOrder>=32]" << endl;
cout << "+ sum (default)" << endl;
cout << "- average" << endl;
return -1;
}

int main(int argC, char* argV[])
{
if (argC<3 || argC>6)
	return Help();

bool sum = true;
bool hasMode = false;
const char* mode = argV[1];
if ((mode[0]=='+' || mode[0]=='-') && mode[1]==0) {
	hasMode = true;
	if (argC<4)
		return Help();
	sum = mode[0]=='+';
	}
else if (argC>5)
	return Help();

const char* iFile = argV[1+hasMode];
const char* oFile = argV[2+hasMode];
	
cout << "Input file: " << iFile << endl;
cout << "Output file: " << oFile << endl;

if (sum) {
	cout << "Summing input values";
	if (!hasMode)
		cout << " (by default)";
	cout << " to the output file";
	cout << endl;
	}
else
	cout << "Writing average values to the output file" << endl;

long nside = 512;
if (argC>=4+hasMode) {
	int n = atoi(argV[3+hasMode]);
	if (n>=1)
		nside = n;
	else
		cout << "Warning: The <nside> parameter should be an integer greater than 0. "<< nside << " assumed" << endl;
	}
	
long rowOrder = 32;
if (argC>=5+hasMode) {
	int n = atoi(argV[4+hasMode]);
	if (n>=1 && n<=nside)
		rowOrder = n;
	else {
		if (n>nside)
			rowOrder = nside;
		cout << "Warning: The <rowOrder> parameter should be an integer between 1 and " << nside << ". "<< rowOrder << " assumed" << endl;
		}
	}
	
AgileMap gm;

if (gm.Read(iFile)==0) {
	cout << "File "  << iFile << " loaded:" << endl;
	cout << "Size: " << gm.Rows() << "x" << gm.Cols() << " = " << gm.Rows()*gm.Cols() << endl;
	cout << "Center l, b: " << gm.GetMapCenterL() << ", " << gm.GetMapCenterB() << endl;
	cout << "Xbin, Ybin: " << gm.GetXbin() << ", " << gm.GetYbin() << endl;
	cout << "X0, Y0: " << gm.GetX0() << ", " << gm.GetY0() << endl;
	cout << "Lonpole: " << gm.GetLonpole() << endl << endl;
	}
else {
	cerr << "Failed opening file cts" << endl;
	return -1;
	}
		
HealpixMap hMap(nside);

int positions = hMap.Count();

int* used = new int[positions];
for (int i=0; i<positions; ++i)
	used[i] = 0;

double l, b;
/// double theta, phi;
long ipring;
for (int row=0; row<gm.Rows(); ++row)
	for (int col=0; col<gm.Cols(); ++col) {
		gm.GetCoords(row, col, &l, &b);
		ipring = hMap.Gal2Hpx(l, b);
		if (ipring>=positions || ipring<0)
			cerr << "ipring(" << row << ", " << col << ")=" << ipring << endl;
		else {
			hMap.AddVal(ipring, gm(row, col));
			++used[ipring];
			}
		}
cerr << "Step 1" << endl;
if (!sum)
	for (int i=0; i<positions; ++i)
		if (used[i]>1 && hMap.Val(i)!=0)
			hMap.Val(i) = hMap.Val(i)/used[i];
cerr << "Step 2" << endl;

int maxUsed = -1;
int minUsed = gm.Rows()*gm.Cols()+1;
int zeros = 0;
int firstUsed = positions+1;
int lastUsed = -1;
for (int i=0; i<positions; ++i) {
	if (maxUsed<used[i])
		maxUsed = used[i];
	if (used[i] && minUsed>used[i])
		minUsed = used[i];
	if (used[i]) {
		if (i>lastUsed)
			lastUsed = i;
		if (i<firstUsed)
			firstUsed = i;
		}
	else
		++zeros;
	}
delete[] used;

cout << "Output map info:" << endl;
cout << "Pixels: " << positions << endl;
cout << "Zeroes: " << zeros << endl;
cout << "Non zeroes: " << positions-zeros << endl;
cout << "Max resolution ratio: " << maxUsed << endl;
cout << "Min resolution ratio: " << minUsed << endl;
cout << "First used pixel: " << firstUsed << endl;
cout << "Last used pixel: " << lastUsed << endl;
cout << "Used pixels range: " << lastUsed-firstUsed+1 << endl << endl;

cout << "Writing file " << oFile << " (nside=" << nside << ", rowOrder=" << rowOrder << ")" << endl;

hMap.Save(oFile, rowOrder);
/// int ret = write_healpix_map(vals, nside, oFile,  0, "G");

HealpixMap hMap2;
if (hMap2.Load(oFile)) {
	cout << "Order=" << hMap2.GetOrder() << endl;
	cout << "NSide=" << hMap2.GetNSide() << endl;
	}
else
	cout << "Not loaded" << endl;

return 0;
}
