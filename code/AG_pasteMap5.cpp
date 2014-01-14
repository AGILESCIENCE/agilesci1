


#include <iostream>
#include <stdlib.h>
#include "SkyMap.h"

using std::cerr;
using std::cout;
using std::endl;



int main(int argC, char* argV[])
{
if (argC<4) {
	cout << "Usage: PasteMap mainMap mapToPaste resultMapName [resultResolution]" << endl;
	return 0;
	}

SkyMap map1;
if (!map1.Load(argV[1])) {
	cerr << "Errors loading file " << argV[1] << endl;
	return -1;
	}
// map1.Print();

SkyMap map2;
if (!map2.Load(argV[2])) {
	cerr << "Errors loading file " << argV[2] << endl;
	return -1;
	}
// map2.Print();

double binSize = -1;
if (argC>4)
	binSize = atof(argV[4]);
SkyMap map3 = PasteMap(map1, map2, binSize);
bool saved = map3.Save(argV[3]);
return saved ? 0 : -1;
}

