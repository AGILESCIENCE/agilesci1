////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG past map
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

