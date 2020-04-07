/* Copyright 2016 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @package TBTKtemp
 *  @file main.cpp
 *  @brief New project
 *
 *  Empty template project.
 *
 *  @author Kristofer Björnson
 */

#include "TBTK/Streams.h"
#include "Calculation.h"
#include <iostream>

using namespace std;
using namespace TBTK;

int main(int argc, char **argv){
	Streams::openLog("Log.txt");

	for(int vz = 1; vz <= 10; vz++)
	{
		string outFile = "vz_" + to_string(vz/10.0) + "_chebychev";
		fstream fileStream;
		fileStream.open(outFile + ".hdf5");
		if(fileStream.fail())
      	{
			Calculation calc(outFile, complex<double>(vz));
			calc.InitModel();
			calc.DoScCalc();
			calc.WriteOutput();
		}
	}



	Streams::closeLog();

	return 0;
}
