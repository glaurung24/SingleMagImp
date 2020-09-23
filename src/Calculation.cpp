#include "Calculation.h"

#include "TBTK/Model.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Property/EigenValues.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/Property/LDOS.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#include "TBTK/FileWriter.h"
#include "TBTK/Streams.h"
#include "TBTK/Array.h"
#include "TBTK/FileReader.h"

#include <complex>
#include<math.h>
#include <cstdlib> 
#include <ctime> 

unsigned int Calculation::system_length;
unsigned int Calculation::system_size;
complex<double> Calculation::mu;
complex<double> Calculation::Vz;
complex<double> Calculation::t;
complex<double> Calculation::delta_start;
complex<double> Calculation::coupling_potential;

Calculation::FunctionDelta Calculation::functionDelta;
Calculation::SelfConsistencyCallback Calculation::selfConsistencyCallback;

const double Calculation::EPS = 1E-4;
const complex<double> Calculation::I = complex<double>(0.0, 1.0);

Array<complex<double>> Calculation::delta;
Array<complex<double>> Calculation::delta_old;

Solver::Diagonalizer Calculation::solver;
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

    
    
    t = 1;
    mu = -0.5; //-1.1, 2.5
    
    delta_start = 0.103229725288; //0.551213123012; //0.0358928467732;
    Vz = vz_input;
    // A coupling potential of 2.5 gives a delta of 0.551213123012
    // A coupling potential of 1.475 gives a delta of 0.103229725288
    coupling_potential = 1.475; //2.0, 1.5 //TODO change back!!!
    delta = Array<complex<double>>({system_size, system_size}, delta_start);

     // //Put random distribution into delta
     srand((unsigned)time(0)); 
     
    for(unsigned int i=0; i<system_size; i++){ 
       for(unsigned int j=0; j<system_size; j++){
        delta[{i, j}] = (complex<double>((rand()%100)+1) + complex<double>((rand()%100)+1)*I)*delta_start/1000.0; 
        
        }
    } 
    delta_old = delta;
    symmetry_on = false;
    use_gpu = false;

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
    

    delta = ConvertVectorToArray(delta_real_from_file, system_size, system_size);
    delta_old = delta;  
    delete [] dims;
    delete [] delta_real_from_file;
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
                        model << HoppingAmplitude(-Vz*2.0*(0.5-s), {x, y, s+2}, {x, y, s+2});
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
    solver.setModel(model);
    solver.setMaxIterations(1000);
    // solver.setScaleFactor(10);
    // solver.setNumCoefficients(1000);
    // solver.setUseLookupTable(true);
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

bool Calculation::SelfConsistencyCallback::selfConsistencyCallback(Solver::Diagonalizer &solver)
{
//    return true; //TODO
    PropertyExtractor::Diagonalizer pe(solver);
    // PropertyExtractor::ChebyshevExpander pe(solver);

    delta_old = delta;
    double diff = 0.0;
    pe.setEnergyWindow(-10, 0, 1000);


    for(unsigned int x=0; x < system_size; x++)
    {
        for(unsigned int y = 0; y < system_size; y++)
        {
            delta[{x , y}] = (-pe.calculateExpectationValue({x,y, 3},{x, y, 0})*coupling_potential*0.5 + delta_old[{x , y}]*0.5);
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
    solver.setModel(model);
    SelfConsistencyCallback selfConsistencyCallback;
    solver.setSelfConsistencyCallback(selfConsistencyCallback);
    solver.setMaxIterations(100);
    solver.run();
	Streams::out << "finished calc" << endl;
}

void Calculation::WriteOutput()
{
	PropertyExtractor::Diagonalizer pe(solver);
  FileWriter::setFileName(outputFileName);

  const double UPPER_BOUND = 4; //10*abs(delta_start);
	const double LOWER_BOUND = -4; //-10*abs(delta_start);
	const int RESOLUTION = 2000;
	pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);



  //Extract DOS and write to file
	Property::DOS dos = pe.calculateDOS();
	FileWriter::writeDOS(dos);

	//Extract eigen values and write these to file
	Property::EigenValues ev = pe.getEigenValues();
	FileWriter::writeEigenValues(ev);

	//Extract LDOS and write to file
	Property::LDOS ldos = pe.calculateLDOS(
		{IDX_X, IDX_Y, IDX_SUM_ALL},
		{system_size, system_size,	4}
	);
	FileWriter::writeLDOS(ldos);

  WriteDelta(0);

  for(unsigned int spin = 0; spin < 4; spin++ )
  {
    string propertyname = "chargeDensity_spin_" + to_string(spin);
    const int dims[2] = {system_size, system_size};
    Array<complex<double>> charge = CalculateChargeDensity(spin);
    FileWriter::write(charge.getData().getData(), 2, dims, propertyname);
  }



//   int nr_excited_states = 30;

//   for(int i = 1; i <= nr_excited_states; i++){
//       Property::WaveFunctions wf = pe.calculateWaveFunctions(
//           {{IDX_ALL, IDX_ALL, IDX_ALL}},
//           {system_size*4+i-nr_excited_states/2}
//       );
//       FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));

//   }

}

void Calculation::WriteDelta(int nr_loop)
{
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