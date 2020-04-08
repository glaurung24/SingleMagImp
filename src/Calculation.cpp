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
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#include "TBTK/PropertyExtractor/ArnoldiIterator.h"
#include "TBTK/Streams.h"
#include "TBTK/Array.h"
#include "TBTK/Exporter.h"
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
double Calculation::energy_bandwidth;
complex<double> Calculation::mu;
complex<double> Calculation::Vz;
complex<double> Calculation::t;
complex<double> Calculation::delta_start;
complex<double> Calculation::coupling_potential;

Calculation::FunctionDelta Calculation::functionDelta;
// Calculation::SelfConsistencyCallback Calculation::selfConsistencyCallback;

const double Calculation::EPS = 1E-4;
const complex<double> Calculation::I = complex<double>(0.0, 1.0);

Array<complex<double>> Calculation::delta;
Array<complex<double>> Calculation::delta_old;

// Solver::Diagonalizer Calculation::solver;
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
    system_length = 10;
    system_size = system_length + 1;
    
    
    t = 1;
    mu = -0.5; //-1.1, 2.5
    
    delta_start = 0.551213123012; //0.0358928467732;
    Vz = vz_input*delta_start/10.0;
    coupling_potential = 2.5; //2.0, 1.5
    delta = Array<complex<double>>({system_size, system_size}, delta_start);

    // // //Put random distribution into delta
    // srand((unsigned)time(0)); 
     
    // for(unsigned int i=0; i<system_size; i++){ 
    //     for(unsigned int j=0; j<system_size; j++){
    //     delta[{i, j}] = (complex<double>((rand()%100)+1) + complex<double>((rand()%100)+1)*I)*delta_start/1000.0; 
         
    //     }
    // } 
    delta_old = delta;
    symmetry_on = false;
    use_gpu = true;
    chebychev_coefficients = 1000;
    energy_points = 2*chebychev_coefficients;
    energy_bandwidth = 8;

    max_arnoldi_iterations = 4000;
    num_eigenvals = 30;
    num_lanczos_vecs = 2*num_eigenvals;

    outputFileName = outputfilename + ".hdf5";
}



void Calculation::InitModel()
{
    //Create model and set up hopping parameters
    model = Model();

    if(!symmetry_on)
    {
        for(unsigned int x = 0; x < system_size; x++){
            for(unsigned int y = 0; y < system_size; y++){
                for(unsigned int s = 0; s < 2; s++){

    //------------------------chemical Potential-----------------------------------
                    //Add hopping amplitudes corresponding to chemical potential
                    model << HoppingAmplitude(-mu,	{x, y, s},	{x, y, s});
                    model << HoppingAmplitude(mu,	{x, y, s+2},	{x, y, s+2});

    //-------------------BCS interaction term------------------------------------------

    //                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
                    model << HoppingAmplitude(Calculation::functionDelta, {x,y,s}, {x,y,(3-s)}) + HC;


    //------------------------Nearest neighbour hopping term--------------------------------------
                    //Add hopping parameters corresponding to t
                    if(x == system_size - 1){
                        model << HoppingAmplitude(-t,	{(x+1)%system_size, y, s},	{x, y, s}) + HC;
                        model << HoppingAmplitude(t,	{x, y, s+2},{(x+1)%system_size, y, s+2}) + HC;
                    }
                    else
                    {
                        model << HoppingAmplitude(-t,	{(x+1)%system_size, y, s},	{x, y, s}) + HC;
                        model << HoppingAmplitude(t,	{(x+1)%system_size, y, s+2},{x, y, s+2}) + HC;
                    }
                    
                    if(y == system_size - 1){
                        model << HoppingAmplitude(-t,	{x, (y+1)%system_size, s},	{x, y, s}) + HC;
                        model << HoppingAmplitude(t,  {x, y, s+2}, {x, (y+1)%system_size, s+2}) + HC;
                    }
                    else
                    {
                        model << HoppingAmplitude(-t,	{x, (y+1)%system_size, s},	{x, y, s}) + HC;
                        model << HoppingAmplitude(t,  {x, y, s+2}, {x, (y+1)%system_size, s+2}) + HC;
                    }
                    

    //---------------------------Zeeman term------------------------------------------
                    if(x == system_size/2 and  y == system_size/2)
                    {
                        model << HoppingAmplitude(Vz*2.0*(0.5-s), {x, y, s}, {x, y, s});
                        model << HoppingAmplitude(Vz*2.0*(0.5-s), {x, y, s+2}, {x, y, s+2});
                    }
                }
            }
        }
    }
    else
    {

    //----------------------using C4 symmetry for calculation--------------------
        int counter = 0;
        for(unsigned int x = 0; x < system_size; x++){
            for(unsigned int y = 0; y < system_size; y++){
                for(unsigned int s = 0; s < 2; s++){

    //------------------------chemical Potential-----------------------------------
                    model << HoppingAmplitude(-mu,	{x, y, s},	{x, y, s});
                    model << HoppingAmplitude(mu,	{x, y, s+2},	{x, y, s+2});

    //-------------------BCS interaction term------------------------------------------

                    model << HoppingAmplitude(Calculation::functionDelta, {x,y,s}, {x,y,(3-s)}) + HC;


    //------------------------Nearest neighbour hopping term--------------------------------------
                    // Add hopping parameters corresponding to t

                    // if(x+1 < system_size){
                    //     if(x == 0 || x+1 == system_size-1){
                    //         model << HoppingAmplitude(-2.0*t,	{x, y, s},	{x+1, y, s}) + HC;
                    //         model << HoppingAmplitude(2.0*t,	{x, y, s+2},{x+1, y, s+2}) + HC;
                    //     }
                    //     // else if(x+1 == system_size){
                    //     //     model << HoppingAmplitude(-t,	{x, y, s},	{x-1, y, s}) + HC;
                    //     //     model << HoppingAmplitude(t,	{x, y, s+2},{x-1, y, s+2}) + HC;
                    //     // }
                    //     else {
                    //         model << HoppingAmplitude(-t,	{x, y, s},	{x+1, y, s}) + HC;
                    //         model << HoppingAmplitude(t,	{x, y, s+2},{x+1, y, s+2}) + HC;
                    //     }
                    // }
                    // if(y+1 < system_size){
                    //     if(y == 0 || (y+1 == system_size-1)){
                    //         model << HoppingAmplitude(-2.0*t,	{x, y, s},	{x, y+1, s}) + HC;
                    //         model << HoppingAmplitude(2.0*t,  {x, y, s+2}, {x, y+1, s+2}) + HC;
                    //     }
                    //     // else if( y+1 == system_size){
                    //     //     model << HoppingAmplitude(-t,	{x, y, s},	{x, y-1, s}) + HC;
                    //     //     model << HoppingAmplitude(t,  {x, y, s+2}, {x, y-1, s+2}) + HC;
                    //     // }
                    //     else {
                    //         model << HoppingAmplitude(-t,	{x, y, s},	{x, y+1, s}) + HC;
                    //         model << HoppingAmplitude(t,  {x, y, s+2}, {x, y+1, s+2}) + HC;
                    //     }
                    // }
                    if(x == 0 || x+1 == system_size){
                        if( x == y)
                        {
                            counter++;
                        }
                        model << HoppingAmplitude(-0.0*t,	{x, y, s},	{x, y, s});
                        model << HoppingAmplitude(0.0*t,	{x, y, s+2},{x, y, s+2});
                    }
                    // else if(x+1 == system_size){
                    //     model << HoppingAmplitude(-t,	{x, y, s},	{x-1, y, s}) + HC;
                    //     model << HoppingAmplitude(t,	{x, y, s+2},{x-1, y, s+2}) + HC;
                    // }
                    if(x+1 < system_size) {
                        model << HoppingAmplitude(-t,	{x, y, s},	{x+1, y, s}) + HC;
                        model << HoppingAmplitude(t,	{x, y, s+2},{x+1, y, s+2}) + HC;
                    }
                    if((y == 0 || y+1 == system_size) && x != 0 && x+1 != system_size){
                        if( x == y)
                        {
                            counter++;
                        }
                        model << HoppingAmplitude(-0.0*t,	{x, y, s},	{x, y, s});
                        model << HoppingAmplitude(0.0*t,  {x, y, s+2}, {x, y, s+2});
                        // cout << "hej x: " << x << ", y: " << y << ", s: " << s << endl;
                    }
                    if(y+1 < system_size) {
                        model << HoppingAmplitude(-t,	{x, y, s},	{x, y+1, s}) + HC;
                        model << HoppingAmplitude(t,  {x, y, s+2}, {x, y+1, s+2}) + HC;
                    }
                    cout << "ups" << endl;
    //---------------------------Zeeman term------------------------------------------
                    // if(x == 0 and  y == 0)
                    // {
                    //     model << HoppingAmplitude(Vz*2.0*(0.5-s), {x, y, s}, {x, y, s});
                    //     model << HoppingAmplitude(Vz*2.0*(0.5-s), {x, y, s+2}, {x, y, s+2});
                    // }
                }
            }
        }
    cout << counter << endl;
    }
    



    model.construct();
}


complex<double> Calculation::FunctionDelta::getHoppingAmplitude(const Index& from, const Index& to) const
{
    unsigned int from_x = from.at(0);
    unsigned int from_y = from.at(1);
    unsigned int from_s = from.at(2);

//    if(isnan(real(delta[from_x])) || isnan(imag(delta[from_x])) ) {
//      cout << "error in " << from_x << endl;
//    }
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

bool Calculation::selfConsistencyCallback(Solver::ChebyshevExpander &solver)
{
    // PropertyExtractor::Diagonalizer pe(solver);
    PropertyExtractor::ChebyshevExpander pe(solver);

	pe.setEnergyWindow(-1*energy_bandwidth, 0, energy_points/2);

    delta_old = delta;
    double diff = 0.0;


    #pragma omp parallel // num_threads( 4 )
    #pragma omp for
    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            delta[{x , y}] = (-pe.calculateExpectationValue({x,y, 3},{x, y, 0})*coupling_potential*0.5 + delta_old[{x , y}]*0.5);
        }
    }


    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            if(abs((delta[{x , y}]-delta_old[{x , y}])/delta_start) > diff)
            {
                diff = abs(delta[{x , y}]-delta_old[{x , y}]);
            }
        }
    }

    diff = diff/abs(delta_start);
    Streams::out << "Updated delta, ddelta = " << to_string(diff) << endl;
    if(diff < EPS)
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
    unsigned int max_iterations = 200;
    solver.setModel(model);
    solver.setScaleFactor(4*energy_bandwidth);
    solver.setNumCoefficients(chebychev_coefficients);
    solver.setUseLookupTable(use_gpu);
    solver.setCalculateCoefficientsOnGPU(use_gpu);
    for(unsigned int loop_counter = 0; loop_counter < max_iterations; loop_counter++)
    {
        if(selfConsistencyCallback(solver))
        {
            break;
        }
    }
	Streams::out << "finished calc" << endl;
}

void Calculation::WriteOutput()
{
	PropertyExtractor::ChebyshevExpander pe(solver);
    Exporter exporter;
    exporter.setFormat(Exporter::Format::ColumnMajor);

    FileWriter::setFileName(outputFileName);

    const double UPPER_BOUND = energy_bandwidth; //10*abs(delta_start);
	const double LOWER_BOUND = -1*energy_bandwidth; //-10*abs(delta_start);
	const int RESOLUTION = energy_points/2;
	pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);

	//Extract LDOS and write to file
	Property::LDOS ldos = pe.calculateLDOS(
		{IDX_X, IDX_Y, IDX_SUM_ALL},
		{system_size, system_size,	4}
	);
    FileWriter::writeLDOS(ldos);

	// exporter.save(ldos, "ldos.csv");

    WriteDelta(0);


    Solver::ArnoldiIterator aSolver;
    aSolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
    aSolver.setModel(model);
    aSolver.setCentralValue(0.0);
    aSolver.setNumEigenValues(num_eigenvals);
    aSolver.setCalculateEigenVectors(true);
    aSolver.setNumLanczosVectors(num_lanczos_vecs);
    aSolver.setMaxIterations(max_arnoldi_iterations);
    aSolver.run();

    PropertyExtractor::ArnoldiIterator ape(aSolver);

	ape.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);

	//Extract eigenvalues and write these to file
	Property::EigenValues ev = ape.getEigenValues();
	FileWriter::writeEigenValues(ev);

	//Extract DOS and write to file
	Property::DOS dos = ape.calculateDOS();
	FileWriter::writeDOS(dos);

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

    // wf = ape->calculateWaveFunction(
	// 	{{IDX_ALL,		IDX_ALL,	IDX_ALL}},
	// 	{NUM_EIGEN_VALUES/2}
	// );
	// FileWriter::writeWaveFunction(wf, "WaveFunctionMFFull");

    // string wave_fct_text = "WaveFunction_";

    // int nr_excited_states = 30;

    // for(int i = 1; i <= nr_excited_states; i++){
    //     wf = ape->calculateWaveFunction(
    //         {{IDX_ALL,		SIZE_Y/2,	IDX_ALL}},
    //         {NUM_EIGEN_VALUES/2+i}
    //     );
    //     FileWriter::writeWaveFunction(wf, "WaveFunction_" + to_string(i));

    //     wf = ape->calculateWaveFunction(
	// 	{{IDX_ALL,		IDX_ALL,	IDX_ALL}},
	// 	{NUM_EIGEN_VALUES/2+i}
    //     );
    //     FileWriter::writeWaveFunction(wf, wave_fct_text + to_string(i) + "Full");
    // }

}

Solver::ArnoldiIterator Calculation::runArnoldiIterator()
{
        Solver::ArnoldiIterator aSolver;
        aSolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
        aSolver.setModel(model);
        aSolver.setCentralValue(0.0);
        aSolver.setNumEigenValues(num_eigenvals);
        aSolver.setCalculateEigenVectors(true);
        aSolver.setNumLanczosVectors(num_lanczos_vecs);
        aSolver.setMaxIterations(max_arnoldi_iterations);
        aSolver.run();
}

void Calculation::WriteDelta(int nr_loop)
{
    // Exporter exporter;
    // exporter.setFormat(Exporter::Format::ColumnMajor);
    // exporter.save(delta, "delta.csv");
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


void Calculation::setOutputFileName(string input)
{
  outputFileName = input;
}

void Calculation::setcoupling_potential(complex<double> input)
{
  coupling_potential = input;
}
