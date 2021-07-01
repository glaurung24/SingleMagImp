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

	// for(int mu = 0; mu <= 8; mu++)
	// {
		// string old_outFile = "";
		string outFile;
		string delta_input_file;
		for(int vz = 0; vz <= 32; vz++)
		{
			// unsigned int nr_phase = 32;
			// for(unsigned int phase_calc = 0; phase_calc <= nr_phase; phase_calc++){
				double Vz = vz/8.0;
				// double Mu = mu/2.0;
				// double phase = static_cast<double>(phase_calc)/nr_phase*M_PI;
				double phase = 0;
				outFile = "vz_" + to_string(Vz) + "mu_" + "-0.5" + "phase_" + to_string(phase) + "_diag_size101_coupling_01_Full_probe30_ldos";
				// outFile = "vz_" + to_string(Vz) + "mu_" + "-0.5" + "_diag_size21noSc";
		        // delta_input_file = "vz_" + to_string(Vz) + "_diag_size21";

				if(!file_exists(outFile + ".hdf5"))
				{
					Calculation calc(outFile, complex<double>(Vz));
					// calc.setMu(Mu);
					// if(old_outFile != "")
					// {
					// calc.readDelta(0, delta_input_file + ".hdf5");
					calc.WriteDelta(0);
					// }
					calc.setPhase(phase);
					calc.InitModel();
					// calc.DoScCalc();
					// calc.WriteOutputSc();
					calc.DoCalc();
					calc.WriteOutput();
					
				}
			// }
			// old_outFile = outFile;
		// }
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
