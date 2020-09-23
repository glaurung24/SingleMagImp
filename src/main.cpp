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

	for(int mu = -8; mu <= 8; mu++)
	{
		string old_outFile = "";
		string outFile;
		for(int vz = 0; vz <= 32; vz++)
		{
			double Vz = vz/16.0;
			double Mu = mu/2.0;
			outFile = "vz_" + to_string(Vz) + "mu_" + to_string(Mu) + "_normalState_diag_size21";
			fstream fileStream;
			fileStream.open(outFile + ".hdf5");
			if(fileStream.fail())
			{
				Calculation calc(outFile, complex<double>(Vz));
				calc.setMu(Mu);
				if(old_outFile != "")
				{
					calc.readDelta(0, old_outFile + ".hdf5");
				}
				calc.InitModel();
				calc.DoScCalc();
				calc.WriteOutput();
			}
			old_outFile = outFile;
		}
	}

	// // for(double coupling = 0.0; coupling <= 3.0; coupling = coupling + 0.05)
	// // {
	// 	string outFile = "coupling_" + to_string(coupling) + "_diag_size21";
	// 	fstream fileStream;
	// 	fileStream.open(outFile + ".hdf5");
	// 	if(fileStream.fail())
    //   	{
	// 		Calculation calc(outFile, coupling);
	// 		calc.InitModel();
	// 		calc.DoScCalc();
	// 		calc.WriteOutput();
	// 	}
	// // }



	Streams::closeLog();

	return 0;
}
