#include "Calculation.h"

#include "TBTK/Model.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/Solver/ArnoldiIterator.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Property/EigenValues.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/Property/LDOS.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
//#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#include "TBTK/PropertyExtractor/ArnoldiIterator.h"
#include "TBTK/Streams.h"
#include "TBTK/Array.h"
#include "TBTK/Exporter.h"
#include "TBTK/FileReader.h"
#include "TBTK/FileWriter.h"

#include <complex>
#include<math.h>
#include <cstdlib> 
#include <ctime> 

unsigned int Calculation::system_length;
unsigned int Calculation::system_size;
unsigned int Calculation::energy_points;
unsigned int Calculation::chebychev_coefficients;
unsigned int Calculation::max_arnoldi_iterations;
unsigned int Calculation::num_eigenvals;
unsigned int Calculation::num_lanczos_vecs;
unsigned int Calculation::max_sc_iterations;
double Calculation::energy_bandwidth;
complex<double> Calculation::mu;
complex<double> Calculation::Vz;
complex<double> Calculation::t;
complex<double> Calculation::delta_start;
complex<double> Calculation::coupling_potential;
complex<double> Calculation::t_probe;
complex<double> Calculation::t_probe_sample;
complex<double> Calculation::t_sample_imp;
double Calculation::phase;
unsigned int Calculation::probe_length;
complex<double> Calculation::delta_probe;
complex<double> Calculation::mu_probe;

Calculation::FunctionDelta Calculation::functionDelta;
Calculation::FunctionDeltaProbe Calculation::functionDeltaProbe;
Calculation::SelfConsistencyCallback Calculation::selfConsistencyCallback;

const double Calculation::EPS = 1E-4;
const complex<double> Calculation::I = complex<double>(0.0, 1.0);

Array<complex<double>> Calculation::delta;
Array<complex<double>> Calculation::delta_old;

Solver::Diagonalizer Calculation::solver;
Solver::ArnoldiIterator Calculation::Asolver;
// Solver::ChebyshevExpander Calculation::solver;

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
    system_length = 30;
    system_size = system_length + 1;

    probe_length = system_length^2;

    delta_start = 0.12188909765277404; // 0.103229725288; //0.551213123012; //0.0358928467732;
    
    t = 1;
    mu = -0.5; //-1.1, 2.5
    t_probe = t;
    t_probe_sample = 0.01*t;
    t_sample_imp = 1.0*t;
    phase = 0;
    delta_probe = delta_start*std::exp(I*phase);
    model_tip = false;
    

    Vz = vz_input;
    // A coupling potential of 2.5 gives a delta of 0.551213123012
    // A coupling potential of 1.475 gives a delta of 0.103229725288
    coupling_potential = 1.475; //2.0, 1.5 //TODO change back!!!
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
    chebychev_coefficients = 1000;
    energy_points = 2*chebychev_coefficients;
    energy_bandwidth = 8;
    max_sc_iterations = 300;

    max_arnoldi_iterations = 4000;
    num_eigenvals = 30;
    num_lanczos_vecs = 2*num_eigenvals;

    outputFileName = outputfilename + ".hdf5";
}

void Calculation::readDelta(int nr_sc_loop, string filename = "")
{
    stringstream loopFileNameReal;

    if(nr_sc_loop < 10)
    {
        loopFileNameReal << "deltaReal" << nr_sc_loop;
    }
    else
    {
        loopFileNameReal << "deltaReal" << nr_sc_loop;
    }
    

    if(filename == "")
    {
        filename = outputFileName;
    }

    
    FileReader::setFileName(filename);
    double* delta_real_from_file = nullptr;
    int rank;
    int* dims;
    FileReader::read(&delta_real_from_file, &rank, &dims, loopFileNameReal.str());
    
    Array<complex<double>> input = ConvertVectorToArray(delta_real_from_file, dims[0], dims[1]);
    delta = deltaPadding(input, system_size, system_size, dims[0], dims[1]);
    delta_old = delta;
    delete [] dims;
    delete [] delta_real_from_file;
}


Array<complex<double>> Calculation::deltaPadding(Array<complex<double>>  input, unsigned int sizeX, unsigned int sizeY, 
                                                unsigned int sizeX_small, unsigned int sizeY_small)
{
    Array<complex<double>> out = Array<complex<double>>({sizeX, sizeX}, delta_start);
    unsigned int marginX = (sizeX - sizeX_small)/2;
    unsigned int marginY = (sizeY - sizeY_small)/2;
    for(unsigned int i=0; i < sizeX_small; i++)
    {
        for(unsigned int j=0; j < sizeY_small; j++)
        {
            out[{i+marginX,j+marginY}] = input[j+i*sizeY_small];
        }
    }
    return out;

}





void Calculation::InitModel()
{
    //Create model and set up hopping parameters
    model = Model();

    // Indeces: sample:0, tip: 1, x, y, spin
    unsigned int system_index_sub = 0; //0 for the substrae, 1 for the tip
    unsigned int system_index_imp = 1;
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
                    model << HoppingAmplitude(Vz*2.0*(0.5-s), {system_index_imp,x, y, s}, {system_index_imp, x, y, s});
                    model << HoppingAmplitude(-Vz*2.0*(0.5-s), {system_index_imp,x, y, s+2}, {system_index_imp,x, y, s+2});
                    model << HoppingAmplitude(-mu,	{system_index_imp,x, y, s},	{system_index_imp,x, y, s});
                    model << HoppingAmplitude(mu,	{system_index_imp,x, y, s+2},	{system_index_imp,x, y, s+2});
                    model << HoppingAmplitude(-t_sample_imp,	{system_index_sub, x, y, s},	{system_index_imp, x, y, s}) + HC;
                    model << HoppingAmplitude(t_sample_imp,  {system_index_sub, x, y, s+2}, {system_index_imp, x, y, s+2}) + HC;
                }
            }
        }
    }
    unsigned int system_index_tip = 2;
    if(model_tip)
    {
        unsigned int position = system_size/2;
        for(unsigned int s = 0; s < 2; s++){
            model << HoppingAmplitude(-t_probe_sample,	{system_index_imp, position, position, s},	{system_index_tip, position, position, s}) + HC;
            model << HoppingAmplitude(t_probe_sample,  {system_index_imp, position, position, s+2}, {system_index_tip, position, position, s+2}) + HC;
        
            for(unsigned pos = system_index_tip; pos < probe_length; pos++){
                if(pos+1 < probe_length){
                    model << HoppingAmplitude(-t_probe,	{pos, position, position, s},	{pos+1, position, position, s}) + HC;
                    model << HoppingAmplitude(t_probe,  {pos, position, position, s+2}, {pos+1, position, position, s+2}) + HC;
                }
                model << HoppingAmplitude(-mu,	{pos, position, position, s},	{pos, position, position, s});
                model << HoppingAmplitude(mu,	{pos, position, position, s+2},	{pos, position, position, s+2});

                // model << HoppingAmplitude(Calculation::functionDeltaProbe, {system_index_tip, pos,s}, {system_index_tip, pos,(3-s)}) + HC;
                model << HoppingAmplitude(Calculation::functionDeltaProbe, {pos, position, position,s}, {pos, position, position,(3-s)}) + HC;
            }
        }
    }
    



    model.construct();
    solver.setModel(model);
    // solver.setMaxIterations(1000);
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
    unsigned int from_s = from.at(3);
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


bool Calculation::SelfConsistencyCallback::selfConsistencyCallback(Solver::Diagonalizer &solver)
{
//    return true; //TODO
    PropertyExtractor::Diagonalizer pe(solver);
    // PropertyExtractor::ChebyshevExpander pe(solver);

	pe.setEnergyWindow(-1*energy_bandwidth, 0, energy_points/2);

    delta_old = delta;
    double diff = 0.0;
    pe.setEnergyWindow(-10, 0, 1000);


    #pragma omp parallel // num_threads( 4 )
    #pragma omp for
    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            delta[{x , y}] = (-pe.calculateExpectationValue({0,x,y, 3},{0,x, y, 0})*coupling_potential*0.7 + delta_old[{x , y}]*0.3);
            if(abs((delta[{x , y}]-delta_old[{x , y}])/delta_start) > diff)
            {
                diff = abs(delta[{x , y}]-delta_old[{x , y}]);
            }
        }
    }

    diff = diff/abs(delta_start);
    Streams::out << "Updated delta, ddelta = " << to_string(diff) << endl;
    if(diff < EPS*abs(delta_start))
    {
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
    // solver.setScaleFactor(4*energy_bandwidth);
    // solver.setNumCoefficients(chebychev_coefficients);
    // solver.setUseLookupTable(use_gpu);
    // solver.setCalculateCoefficientsOnGPU(use_gpu);
    cout << "@sc calc" << system_size  << endl;
    for(unsigned int loop_counter = 0; loop_counter < max_sc_iterations; loop_counter++)
    {
        if(selfConsistencyCallback.selfConsistencyCallback(solver))
        {
            break;
        }
    }
	Streams::out << "finished calc" << endl;
}

void Calculation::DoCalc()
{
    model.construct();
    // Asolver.setModel(model);
    // Asolver.setNumLanczosVectors(4000);
    // Asolver.setMaxIterations(20000);
    // Asolver.setNumEigenValues(2000);
    // Asolver.setCalculateEigenVectors(false);
    // Asolver.setCentralValue(-0.01);
    // Asolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
    // Asolver.run();
    solver.setModel(model);
    solver.run();
	Streams::out << "finished calc" << endl;
}

void Calculation::WriteOutputSc()
{
    // PropertyExtractor::ChebyshevExpander pe(solver);
    // Exporter exporter;
    // exporter.setFormat(Exporter::Format::ColumnMajor);
    PropertyExtractor::Diagonalizer pe(solver);

    // const double UPPER_BOUND = 5; //10*abs(delta_start);
	// const double LOWER_BOUND = -5; //-10*abs(delta_start);
	// const int RESOLUTION = 2000;
	// pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);


    FileWriter::setFileName(outputFileName);

  //Extract DOS and write to file
	Property::DOS dos = pe.calculateDOS();
	FileWriter::writeDOS(dos);

    cout << "@write out" << system_size  << endl;

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

	// exporter.save(ldos, "ldos.csv");

  WriteDelta(0);

//   for(unsigned int spin = 0; spin < 4; spin++ )
//   {
//     string propertyname = "chargeDensity_spin_" + to_string(spin);
//     const int dims[2] = {system_size, system_size};
//     Array<complex<double>> charge = CalculateChargeDensity(spin);
//     FileWriter::write(charge.getData().getData(), 2, dims, propertyname);
//   }



  int nr_excited_states = 30;

  for(int i = 1; i <= nr_excited_states; i++){
      Property::WaveFunctions wf = pe.calculateWaveFunctions(
          {{0, IDX_ALL, IDX_ALL, IDX_ALL}},
          {system_size*4+i-nr_excited_states/2}
      );
      FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));
   }
}


void Calculation::WriteOutput()
{
    PropertyExtractor::ArnoldiIterator pe(Asolver);
    // PropertyExtractor::Diagonalizer pe(solver);
    FileWriter::setFileName(outputFileName);

    // const double UPPER_BOUND = 2*abs(delta_start);
	// const double LOWER_BOUND = -2*abs(delta_start);
	// const int RESOLUTION = 5000;
	// pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);

// //Extract spin-polarized LDOS and write to file
// 	Property::SpinPolarizedLDOS spLdos = ape->calculateSpinPolarizedLDOS(
// 		{IDX_X,		IDX_Y,	IDX_SPIN},
// 		{system_size, system_size, 4}
// 	);
// 	FileWriter::writeSpinPolarizedLDOS(spLdos);

	// Property::WaveFunction wf = ape->calculateWaveFunction(
	// 	{{IDX_ALL,		SIZE_Y/2,	IDX_ALL}},
	// 	{NUM_EIGEN_VALUES/2}
	// );
	// FileWriter::writeWaveFunction(wf, "WaveFunctionMF");

  //Extract DOS and write to file
	// Property::DOS dos = pe.calculateDOS();
	// FileWriter::writeDOS(dos);

	//Extract eigen values and write these to file
	Property::EigenValues ev = pe.getEigenValues();
	FileWriter::writeEigenValues(ev);

	// Extract LDOS and write to file

        Property::LDOS ldos = pe.calculateLDOS(
            {0, IDX_X, IDX_Y, IDX_SUM_ALL},
            {0, system_size, system_size,	4}
        );
        FileWriter::writeLDOS(ldos);


  int nr_excited_states = num_eigenvals;

  for(int i = 1; i <= nr_excited_states; i++){
      Property::WaveFunctions wf = pe.calculateWaveFunctions(
          {{0, IDX_ALL, IDX_ALL, IDX_ALL}},
          {system_size*4+i-nr_excited_states/2}
      );
      FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));

  }


//   WriteDelta(0);
}

void Calculation::runArnoldiIterator()
{
        Asolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
        Asolver.setModel(model);
        Asolver.setCentralValue(0.0);
        Asolver.setNumEigenValues(num_eigenvals);
        Asolver.setCalculateEigenVectors(true);
        Asolver.setNumLanczosVectors(num_lanczos_vecs);
        Asolver.setMaxIterations(max_arnoldi_iterations);
        Asolver.run();
}

void Calculation::WriteDelta(int nr_loop)
{
    // Exporter exporter;
    // exporter.setFormat(Exporter::Format::ColumnMajor);
    // exporter.save(delta, "delta.csv");
    cout << "@write out delta" << system_size  << endl;
    FileWriter::setFileName(outputFileName);
    const int RANK = 2;
    int dims[RANK] = {system_size, system_size};
    FileWriter::write(GetRealVec(delta).getData().getData(), RANK, dims, "deltaReal" + to_string(nr_loop));
    FileWriter::write(GetImagVec(delta).getData().getData(), RANK, dims, "deltaImag" + to_string(nr_loop));
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

Array<complex<double>> Calculation::CalculateChargeDensity(unsigned int spin)
{
    PropertyExtractor::Diagonalizer pe(solver);
    pe.setEnergyWindow(-10, 0, 1000);
    Array<complex<double>> charge({system_size, system_size});
    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            charge[{x , y}] = pe.calculateExpectationValue({x,y, spin},{x, y, spin});
        }
    }
    return charge;
}


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


void Calculation::setTImpSample(complex<double> input)
{
  t_sample_imp = input;
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
