#include "Calculation.h"

#include "TBTK/Model.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/ChebyshevExpander.h"
// #include "TBTK/Solver/ArnoldiIterator.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Property/EigenValues.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/Property/LDOS.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
//#include "TBTK/PropertyExtractor/ArnoldiIterator.h"
#include "TBTK/Streams.h"
#include "TBTK/Array.h"
#include "TBTK/Exporter.h"
#include "TBTK/Resource.h"

#include <complex>
#include<math.h>
#include <cstdlib> 
#include <ctime> 

unsigned int Calculation::system_length;
unsigned int Calculation::system_size;
unsigned int Calculation::num_chebyshev_coeff;
unsigned int Calculation::num_energy_points;
double Calculation::lower_energy_bound;
double Calculation::upper_energy_bound;
complex<double> Calculation::mu;
complex<double> Calculation::Vz;
complex<double> Calculation::t;
complex<double> Calculation::delta_start;
complex<double> Calculation::coupling_potential;
complex<double> Calculation::t_probe;
complex<double> Calculation::t_probe_sample;
double Calculation::phase;
unsigned int Calculation::probe_length;
complex<double> Calculation::delta_probe;
complex<double> Calculation::mu_probe;

Calculation::FunctionDelta Calculation::functionDelta;
Calculation::FunctionDeltaProbe Calculation::functionDeltaProbe;
// Calculation::SelfConsistencyCallback Calculation::selfConsistencyCallback;

const double Calculation::EPS = 1E-4;
const complex<double> Calculation::I = complex<double>(0.0, 1.0);

Array<complex<double>> Calculation::delta;
Array<complex<double>> Calculation::delta_old;

// Solver::Diagonalizer Calculation::solver;
// Solver::ArnoldiIterator Calculation::Asolver;
Solver::ChebyshevExpander Calculation::solver;

string Calculation::outputFileName = "";

Calculation::Calculation()
{
    //ctor
    Init(string("out") );
}

Calculation::Calculation(string outputfilename, complex<double> vz_input)
{
    //ctor
    Init(outputfilename, vz_input);
}

Calculation::~Calculation()
{
    //dtor
}

void Calculation::Init(string outputfilename, complex<double> vz_input)
{
    system_length = 250;
    system_size = system_length + 1;

    probe_length = system_size^2;

    delta_start = 0.0379431;// 0.12188909765277404; // 0.103229725288; //0.551213123012; //0.0358928467732;
    
    t = 1;
    mu = -0.5; //-1.1, 2.5
    t_probe = t;
    t_probe_sample = 0.01*t;
    phase = 0;
    delta_probe = delta_start*std::exp(I*phase);
    model_tip = false;
    

    Vz = vz_input;
    // A coupling potential of 2.5 gives a delta of 0.551213123012
    // A coupling potential of 1.475 gives a delta of 0.103229725288
    coupling_potential = 1.2;//1.475; //2.0, 1.5 //TODO change back!!!
    delta = Array<complex<double>>({system_size, system_size}, delta_start);

     // //Put random distribution into delta
    //  srand((unsigned)time(0)); 
     
    // for(unsigned int i=0; i<system_size; i++){ 
    //    for(unsigned int j=0; j<system_size; j++){
    //     delta[{i, j}] = (complex<double>((rand()%100)+1) + complex<double>((rand()%100)+1)*I)*delta_start; 
        
    //     }
    // } 

    delta_old = delta;
    symmetry_on = false;
    use_gpu = false;
    num_chebyshev_coeff = 20000;
    num_energy_points = num_chebyshev_coeff * 2;
    lower_energy_bound = -6;
    upper_energy_bound = 6;

    outputFileName = outputfilename;
}

void Calculation::setSystem_length(unsigned int lenght)
{
    system_length = lenght;
    system_size = system_length + 1;
    delta = Array<complex<double>>({system_size, system_size}, delta_start);
}

void Calculation::readDelta(int nr_sc_loop, string filename = "")
{
    if( filename == "")
    {
        filename = DeltaOutputFilename(nr_sc_loop) + ".json";
    }

    Resource resource;
    resource.read(filename);
    delta = Array<complex<double>>(resource.getData(), Serializable::Mode::JSON);
    delta_old = delta;
    unsigned int position = system_size/2;
    delta_start = delta[{position,position}];
}


string Calculation::DeltaOutputFilename(const int nr_sc_loop)
{
    string nr_padded = to_string(nr_sc_loop);
    nr_padded.insert(0, 3 - nr_padded.length(), '0');
    return outputFileName + "_delta_" + nr_padded;
}





void Calculation::InitModel()
{
    //Create model and set up hopping parameters
    model = Model();

    // Indeces: sample:0, tip: 1, x, y, spin
    unsigned int system_index_sub = 0; //0 for the substrae, 1 for the tip
    for(unsigned int x = 0; x < system_size; x++){
        for(unsigned int y = 0; y < system_size; y++){
            for(unsigned int s = 0; s < 2; s++){

//------------------------chemical Potential-----------------------------------
                //Add hopping amplitudes corresponding to chemical potential
                model << HoppingAmplitude(-mu,	{system_index_sub,x, y, s},	{system_index_sub,x, y, s});
                model << HoppingAmplitude(mu,	{system_index_sub,x, y, s+2},	{system_index_sub,x, y, s+2});

//-------------------BCS interaction term------------------------------------------

//                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
                model << HoppingAmplitude(Calculation::functionDelta, {system_index_sub,x,y,s}, {system_index_sub,x,y,(3-s)}) + HC;


//------------------------Nearest neighbour hopping term--------------------------------------
                //Add hopping parameters corresponding to t
                if(x == system_size - 1){
                    model << HoppingAmplitude(-t,	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
                    model << HoppingAmplitude(t,	{system_index_sub,x, y, s+2},{system_index_sub,(x+1)%system_size, y, s+2}) + HC;
                }
                else
                {
                    model << HoppingAmplitude(-t,	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
                    model << HoppingAmplitude(t,	{system_index_sub,(x+1)%system_size, y, s+2},{system_index_sub,x, y, s+2}) + HC;
                }
                
                if(y == system_size - 1){
                    model << HoppingAmplitude(-t,	{system_index_sub, x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
                    model << HoppingAmplitude(t,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
                }
                else
                {
                    model << HoppingAmplitude(-t,	{system_index_sub,x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
                    model << HoppingAmplitude(t,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
                }
                

//---------------------------Zeeman term------------------------------------------
                if(x == system_size/2 and  y == system_size/2)
                {
                    model << HoppingAmplitude(Vz*2.0*(0.5-s), {system_index_sub,x, y, s}, {system_index_sub, x, y, s});
                    model << HoppingAmplitude(-Vz*2.0*(0.5-s), {system_index_sub,x, y, s+2}, {system_index_sub,x, y, s+2});
                }
            }
        }
    }
    unsigned int system_index_tip = 1;
    if(model_tip)
    {
        unsigned int position = system_size/2;
        for(unsigned int s = 0; s < 2; s++){
            model << HoppingAmplitude(-t_probe_sample,	{system_index_sub, position, position, s},	{system_index_tip, 0, s}) + HC;
            model << HoppingAmplitude(t_probe_sample,  {system_index_sub, position, position, s+2}, {system_index_tip, 0, s+2}) + HC;
        
            for(unsigned pos = 0; pos < probe_length; pos++){
                if(pos+1 < probe_length){
                    model << HoppingAmplitude(-t_probe,	{system_index_tip, pos, s},	{system_index_tip, pos+1, s}) + HC;
                    model << HoppingAmplitude(t_probe,  {system_index_tip, pos, s+2}, {system_index_tip, pos+1, s+2}) + HC;
                }
                model << HoppingAmplitude(-mu,	{system_index_tip, pos, s},	{system_index_tip, pos, s});
                model << HoppingAmplitude(mu,	{system_index_tip, pos, s+2},	{system_index_tip, pos, s+2});

                model << HoppingAmplitude(Calculation::functionDeltaProbe, {system_index_tip, pos,s}, {system_index_tip, pos,(3-s)}) + HC;
            }
        }
    }
    



    // solver.setScaleFactor(10);
    // solver.setNumCoefficients(1000);
    // solver.setUseLookupTable(true);
}


complex<double> Calculation::FunctionDelta::getHoppingAmplitude(const Index& from, const Index& to) const
{
    unsigned int from_x = from.at(1);
    unsigned int from_y = from.at(2);
    unsigned int from_s = from.at(3);

    switch(from_s)
    {
    case 0:
        return conj(delta[{from_x, from_y}]);
    case 1:
        return -conj(delta[{from_x, from_y}]);
    case 2:
        return -delta[{from_x, from_y}];
    case 3:
        return delta[{from_x, from_y}];
    default:
        Streams::err << "something went wrong in Calculation::FuncDelta." << endl;
        return 0;
    }
}

complex<double> Calculation::FunctionDeltaProbe::getHoppingAmplitude(const Index& from, const Index& to) const
{
    unsigned int from_s = from.at(2);
    switch(from_s)
    {
    case 0:
        return conj(delta_probe);
    case 1:
        return -conj(delta_probe);
    case 2:
        return -delta_probe;
    case 3:
        return delta_probe;
    default:
        Streams::err << "something went wrong in Calculation::FuncDelta." << endl;
        return 0;
    }
}


// bool Calculation::SelfConsistencyCallback::selfConsistencyCallback(Solver::Diagonalizer &solver)
bool Calculation::selfConsistencyCallback()
{
    // PropertyExtractor::Diagonalizer pe(solver);
    // solver.setVerbose(true); 
    PropertyExtractor::ChebyshevExpander pe(solver);

    delta_old = delta;
    // Array<complex<double>> delta_temp = delta;
    double diff = 0.0;
    pe.setEnergyWindow(lower_energy_bound, 0, num_energy_points/2);


    // for(unsigned int x=0; x < system_size; x++)
    // {
    //     #pragma omp parallel for
    //     for(unsigned int y = 0; y < system_size; y++)
    //     {
    //         delta_temp[{x , y}] = (-pe.calculateExpectationValue({0,x,y, 3},{0,x, y, 0})*coupling_potential*0.5 + delta_old[{x , y}]*0.5);
    //         if(abs((delta_temp[{x , y}]-delta_old[{x , y}])/delta_start) > diff)
    //         {
    //             diff = abs(delta_temp[{x , y}]-delta_old[{x , y}]);
    //         }
    //     }
    // }
    // delta = delta_temp;
    complex<double> delta_temp;
    unsigned int position = system_size/2; 
    delta_temp = -pe.calculateExpectationValue({0,position,position, 3},{0,position, position, 0})*coupling_potential;
    diff = abs(delta_temp-delta_old[{position , position}]);
    diff = diff/abs(delta_start);
    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            delta[{x , y}] = delta_temp;
        }
    }


    Streams::out << "Updated delta = " << to_string(real(delta_temp)) << ", ddelta = " << to_string(diff) << endl;
    if(diff < EPS)
    {
        cout << "finished self consistency loop" << endl;
        return true;
    }
    else
    {
        return false;
    }
}

void Calculation::DoScCalc()
{
    model.construct();
    solver.setModel(model);
    solver.setScaleFactor(upper_energy_bound);
    solver.setCalculateCoefficientsOnGPU(use_gpu);
    solver.setGenerateGreensFunctionsOnGPU(false);
    solver.setUseLookupTable(true);
    solver.setNumCoefficients(num_chebyshev_coeff);
    for(unsigned int i = 0; i < 200; i++)
    {
        if(selfConsistencyCallback())
        {
            break;
        }
    }
    // solver.setSelfConsistencyCallback(selfConsistencyCallback);
    // solver.setMaxIterations(100);
    // solver.run();
	Streams::out << "finished calc" << endl;
}

// void Calculation::DoCalc()
// {
//     model.construct();
//     Asolver.setModel(model);
//     Asolver.setNumLanczosVectors(1600);
//     Asolver.setMaxIterations(5000);
//     Asolver.setNumEigenValues(800);
//     Asolver.setCalculateEigenVectors(true);
//     Asolver.setCentralValue(-0.01);
//     Asolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
//     Asolver.run();
// 	Streams::out << "finished calc" << endl;
// }

void Calculation::WriteOutputSc()
{
	// PropertyExtractor::Diagonalizer pe(solver);
    // PropertyExtractor::ChebyshevExpander pe(solver);
	// pe.setEnergyWindow(lower_energy_bound, upper_energy_bound, num_energy_points);
    WriteDelta(0);
    // Exporter exporter;

  //Extract DOS and write to file
	// Property::DOS dos = pe.calculateDOS();
	// exporter.save(dos, outputFileName + "_dos.csv");

	//Extract eigen values and write these to file
	// Property::EigenValues ev = pe.getEigenValues();
	// exporter.save(ev, outputFileName + "_eigenvalues.csv");

	// Extract LDOS and write to file
    // if(model_tip){
    //     Property::LDOS ldos = pe.calculateLDOS(
    //         {IDX_Z,IDX_X, IDX_Y, IDX_SUM_ALL},
    //         {probe_length,system_size, system_size,	4}
    //     );
    //     FileWriter::writeLDOS(ldos);
    // }
    // else{
        // Property::LDOS ldos = pe.calculateLDOS(
        //     {IDX_Z, IDX_X, IDX_Y, IDX_SUM_ALL},
        //     {1, system_size, system_size,	4}
        // );
        // FileWriter::writeLDOS(ldos);
    // }


//   WriteDelta(0);

//   for(unsigned int spin = 0; spin < 4; spin++ )
//   {
//     string propertyname = "chargeDensity_spin_" + to_string(spin);
//     const int dims[2] = {system_size, system_size};
//     Array<complex<double>> charge = CalculateChargeDensity(spin);
//     FileWriter::write(charge.getData().getData(), 2, dims, propertyname);
//   }



//   int nr_excited_states = 30;

//   for(int i = 1; i <= nr_excited_states; i++){
//       Property::WaveFunctions wf = pe.calculateWaveFunctions(
//           {{IDX_ALL, IDX_ALL, IDX_ALL}},
//           {system_size*4+i-nr_excited_states/2}
//       );
//       FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));

//   }

}

// void Calculation::WriteOutput()
// {
// 	PropertyExtractor::ArnoldiIterator pe(Asolver);
//     FileWriter::setFileName(outputFileName);

//     const double UPPER_BOUND = 1; //10*abs(delta_start);
// 	const double LOWER_BOUND = -1; //10*abs(delta_start);
// 	const int RESOLUTION = 2000;
// 	pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);



//   //Extract DOS and write to file
// 	Property::DOS dos = pe.calculateDOS();
// 	FileWriter::writeDOS(dos);

// 	//Extract eigen values and write these to file
// 	Property::EigenValues ev = pe.getEigenValues();
// 	FileWriter::writeEigenValues(ev);

// 	// Extract LDOS and write to file

//         // Property::LDOS ldos = pe.calculateLDOS(
//         //     {IDX_Z, IDX_X, IDX_Y, IDX_SUM_ALL},
//         //     {1, system_size, system_size,	4}
//         // );
//         // FileWriter::writeLDOS(ldos);



//   WriteDelta(0);
// }

void Calculation::WriteDelta(int nr_loop)
{
    string filename = DeltaOutputFilename(nr_loop);
    Exporter exporter;
    exporter.save(GetRealVec(delta), filename + ".csv" );

    if(nr_loop == 0)
    {
        string delta_out = delta.serialize(Serializable::Mode::JSON);
        Resource resource;
        resource.setData(delta_out);
        resource.write(filename + ".json");
    }


    // FileWriter::setFileName(outputFileName);
    // const int RANK = 2;
    // int dims[RANK] = {system_size, system_size};
    // FileWriter::write(GetRealVec(delta).getData().getData(), RANK, dims, "deltaReal" + to_string(nr_loop));
    // FileWriter::write(GetImagVec(delta).getData().getData(), RANK, dims, "deltaImag" + to_string(nr_loop));

}

Array<double> Calculation::GetRealVec(Array<complex<double>> input)
{
    unsigned size_x = input.getRanges()[0];
    unsigned size_y = input.getRanges()[1];
    Array<double> output({size_x, size_y});

    for(unsigned int i=0; i < size_x; i++)
    {
        for(unsigned int j=0; j < size_y; j++)
        {
            output[{i,j}] = real(input[{i,j}]);
        }
    }
    return output;
}

// Array<complex<double>> Calculation::CalculateChargeDensity(unsigned int spin)
// {
//     PropertyExtractor::Diagonalizer pe(solver);
//     pe.setEnergyWindow(-10, 0, 1000);
//     Array<complex<double>> charge({system_size, system_size});
//     for(unsigned int x=0; x < system_size; x++)
//     {
//         for(unsigned int y = 0; y < system_size; y++)
//         {
//             charge[{x , y}] = pe.calculateExpectationValue({x,y, spin},{x, y, spin});
//         }
//     }
//     return charge;
// }


Array<double> Calculation::GetImagVec(Array<complex<double>> input)
{
    unsigned size_x = input.getRanges()[0];
    unsigned size_y = input.getRanges()[1];
    Array<double> output({size_x, size_y});

    for(unsigned int i=0; i < size_x; i++)
    {
        for(unsigned int j=0; j < size_y; j++)
        {
            output[{i,j}] = imag(input[{i,j}]);
        }
    }
    return output;
}




void Calculation::setVz(complex<double> input)
{
  Vz = input;
}

void Calculation::setMu(complex<double> input)
{
  mu = input;
}

void Calculation::setPhase(double input)
{
  phase = input;
  delta_probe = delta_start*std::exp(I*phase);
}


void Calculation::setOutputFileName(string input)
{
  outputFileName = input;
}

void Calculation::setcoupling_potential(complex<double> input)
{
  coupling_potential = input;
}

Array<complex<double>> Calculation::ConvertVectorToArray(const double *input, unsigned int sizeX, unsigned int sizeY)
{
    Array<complex<double>> out = Array<complex<double>>({sizeX, sizeX}, 0);
    
    for(unsigned int i=0; i < sizeX; i++)
    {
        for(unsigned int j=0; j < sizeY; j++)
        {
            out[{i,j}] = input[j+i*sizeY];
        }
    }
    return out;
}