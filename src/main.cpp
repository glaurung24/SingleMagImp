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
#include <string>
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
		// string outFile2;
		// string delta_input_file;
		// for(int pos = 0; pos <= 5; pos++)
		// {
		for(int vz = 0; vz <= 128; vz++)
		{
			// unsigned int nr_phase = 64;
			//   for(unsigned int phase_calc = 0; phase_calc <= nr_phase; phase_calc++){
			// 		double phase = static_cast<double>(phase_calc)/nr_phase*M_PI;
					// double phase = 0;
					// int pos = 0;
					double Vz = vz/64.0;
					// double Vz = 1.599609;//vz/32.0;

					// outFile = "vz_" + to_string(Vz) +  "mu_" + "-0.5"  + "phase_" + to_string(phase) + "_diag_size151_probeNewPos_" + to_string(pos); 
					outFile = "vz_" + to_string(Vz) +  "Eigenstates_calc"; 
					
					// outFile = "vz_" + to_string(Vz) + "mu_" + "-0.5" + "_diag_size51_sc";
					string delta_input_file =  "vz_" + to_string(Vz) +  "mu_" + "-0.5" + "phase_0.000000" + "_diag_size151_probeNewPos_0_" + "delta_000.csv";
					if(!file_exists(outFile) & file_exists(delta_input_file))
					{
						ofstream output(outFile);
						Calculation calc(outFile, complex<double>(Vz));
						// calc.setTipPosition(pos);
						// calc.setMu(Mu);
						// if(old_outFile != "")
						// {
						calc.readDeltaCsv(0, delta_input_file);
						// unsigned int position = calc.getSystemSize()/2;
						// calc.setTipPosition(position);
						// }
						// calc.setPhase(phase);
						// complex<double> delta_start = calc.getDeltaStart();
						// complex<double> dD = delta_start*complex<double>(ddelta)/10.0;
						// calc.setDeltaDelta(dD);
						//  cout << to_string(real(dD)) << endl;
						calc.InitModel();
						// calc.DoScCalc();
						// calc.WriteOutputSc();
						// calc.DoCalc();
						// calc.WriteOutput();
						calc.CalcEigenstates();
						// calc.WriteDelta(0);
						
					}
			//}
			// old_outFile = outFile;
		}
		// }
	// }
	Streams::closeLog();

	return 0;
}

