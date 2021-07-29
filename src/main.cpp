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
#include <math.h>
#include <sys/stat.h>

using namespace std;
using namespace TBTK;

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char **argv){
	Streams::openLog("Log.txt");


	
	// string outfile_old = "";
	// for(double coupling = 1.2; coupling <= 1.3; coupling = coupling+0.01)
	// {
	// 	unsigned int size = 250;
	// 	string outFile = "coupling_" + to_string(coupling) + "_cheb_size" + to_string(size + 1);
	// 	Calculation calc(outFile, 0);
		
	// 	fstream fileStream;
	// 	fileStream.open( calc.DeltaOutputFilename(0) + ".json");
	// 	if(fileStream.fail())
    //   	{
	// 		fstream fileStreamOld;
	// 		fileStreamOld.open( outfile_old);
	// 		if(!fileStreamOld.fail())
	// 		{
	// 			calc.readDelta(0, outfile_old);
	// 		}
	// 		calc.setcoupling_potential(coupling);
	// 		calc.InitModel();
	// 		calc.DoScCalc();
	// 		calc.WriteOutputSc();
	// 	}
	// 	outfile_old = calc.DeltaOutputFilename(0) + ".json";
	// }


	string outFile;
	string outFile_old = "";
	string delta_input_file;
	for(int vz = 0; vz <= 48; vz++)
	{

		double Vz = vz/16.0;
		outFile = "vz_" + to_string(Vz) + "mu_" + "-0.5"  + "_cheby_size251_sc";
		delta_input_file = "vz_" + to_string(Vz) + "_diag_size21";

		if(!file_exists(outFile + "_delta_000.json"))
		{
			Calculation calc(outFile, complex<double>(Vz));
			cout << outFile_old << endl;
			if(file_exists(outFile_old))
			{
				calc.readDelta(0, outFile_old);
			}
			calc.InitModel();
			calc.DoScCalc();
			calc.WriteOutputSc();
			outFile_old = calc.DeltaOutputFilename(0) + ".json";
		}
		
	}




	Streams::closeLog();

	return 0;
}

