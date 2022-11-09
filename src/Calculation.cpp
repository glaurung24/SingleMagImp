#include "Calculation.h"

#include "TBTK/Model.h"
#ifdef GPU_CALCULATION
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#else
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#endif

#include "TBTK/Solver/ArnoldiIterator.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/Property/EigenValues.h"
#include "TBTK/Property/WaveFunctions.h"
#include "TBTK/Property/LDOS.h"
#include "TBTK/PropertyExtractor/ArnoldiIterator.h"
#include "TBTK/Streams.h"
#include "TBTK/Array.h"
#include "TBTK/Exporter.h"
#include "TBTK/Smooth.h"
// #include "TBTK/Resource.h"
// #include "TBTK/FileReader.h"
// #include "TBTK/FileWriter.h"


#include <complex>
#include <math.h>
#include <cstdlib> 
#include <ctime> 
#include <vector>
#include <fstream>

unsigned int Calculation::system_length;
unsigned int Calculation::system_size;
unsigned int Calculation::delta_simulation_size;
unsigned int Calculation::energy_points;
unsigned int Calculation::chebychev_coefficients;
unsigned int Calculation::max_arnoldi_iterations;
unsigned int Calculation::num_eigenvals;
unsigned int Calculation::num_lanczos_vecs;
unsigned int Calculation::max_sc_iterations;
double Calculation::energy_bandwidth;
unsigned int Calculation::tip_position;
complex<double> Calculation::mu;
complex<double> Calculation::Vz;
complex<double> Calculation::t;
complex<double> Calculation::alpha;
complex<double> Calculation::delta_p;
complex<double> Calculation::delta_s_bond;
complex<double> Calculation::delta_start;
complex<double> Calculation::delta_Delta;
complex<double> Calculation::coupling_potential;
complex<double> Calculation::t_probe;
complex<double> Calculation::t_probe_sample;
complex<double> Calculation::t_sample_imp;
double Calculation::phase;
unsigned int Calculation::probe_length;
complex<double> Calculation::delta_probe;
complex<double> Calculation::mu_probe;

Calculation::FunctionDelta Calculation::functionDelta;
// Calculation::FunctionDeltaProbe Calculation::functionDeltaProbe;
// Calculation::FunctionDeltaDelta Calculation::functionDeltaDelta;
Calculation::SelfConsistencyCallback Calculation::selfConsistencyCallback;

const double Calculation::EPS = 1E-4;
const complex<double> Calculation::I = complex<double>(0.0, 1.0);

TBTK::Array<complex<double>> Calculation::delta;
Array<complex<double>> Calculation::delta_old;

bool Calculation::symmetry_on;
bool Calculation::use_gpu;
bool Calculation::model_tip;
bool Calculation::flat_tip;
bool Calculation::p_wave_sc;

Solver::ArnoldiIterator Calculation::Asolver;

#ifdef GPU_CALCULATION
Solver::ChebyshevExpander Calculation::solver;
#else
Solver::Diagonalizer Calculation::solver;
#endif

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
    system_length = 100;
    delta_simulation_size = 101;
    system_size = system_length + 1;

    probe_length = system_size^2;

    delta_start = 0.115490; //0.12188909765277404; // 0.103229725288; //0.551213123012; //0.0358928467732;
    delta_p = 0.1*delta_start;
    delta_s_bond = 0.1*delta_start;
    
    t = 1;
    mu = -0.5; //-1.1, 2.5
    alpha = 0.0;
    t_probe = t;
    t_probe_sample = 0.1*t;
    t_sample_imp = 1.0*t;
    phase = 0;
    delta_Delta = 0;
    delta_probe = delta_start*std::exp(I*phase);
    model_tip = false;
    flat_tip = false;
    model_hubbard_model = false;
    calculate_waveFcts = true;
    p_wave_sc = true;

    tip_position = system_size/2;
    

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
    #ifdef GPU_CALCULATION
    use_gpu = true;
    chebychev_coefficients = 20000; //old: 25000
    energy_points = chebychev_coefficients * 2;
    #else
    use_gpu = false;
    #endif

    energy_bandwidth = 8;
    max_sc_iterations = 300;

    max_arnoldi_iterations = 20000;
    num_eigenvals = 1500;
    num_lanczos_vecs = 2*num_eigenvals;

    outputFileName = outputfilename;
}


void Calculation::setSystem_length(unsigned int lenght)
{
    system_length = lenght;
    system_size = system_length + 1;
    delta = Array<complex<double>>({system_size, system_size}, delta_start);
}

// void Calculation::readDeltaHdf5(int nr_sc_loop, string filename = "")
// {
//     stringstream loopFileNameReal;

//     if(nr_sc_loop < 10)
//     {
//         loopFileNameReal << "deltaReal" << nr_sc_loop;
//     }
//     else
//     {
//         loopFileNameReal << "deltaReal" << nr_sc_loop;
//     }
    

//     if(filename == "")
//     {
//         filename = outputFileName;
//     }

    
//     FileReader::setFileName(filename + ".hdf5");
//     double* delta_real_from_file = nullptr;
//     int rank;
//     int* dims;
//     FileReader::read(&delta_real_from_file, &rank, &dims, loopFileNameReal.str());
    
//     Array<complex<double>> input = ConvertVectorToArray(delta_real_from_file, dims[0], dims[1]);
//     delta = deltaPadding(input, system_size, system_size, dims[0], dims[1]);
//     delta_old = delta;
//     delete [] dims;
//     delete [] delta_real_from_file;
// }

// void Calculation::readDeltaJson(int nr_sc_loop, string filename = "")
// {
//     if( filename == "")
//     {
//         filename = DeltaOutputFilename(nr_sc_loop, filename) + ".json";
//     }

//     Resource resource;
//     resource.read(filename);
//     cout << resource.getData() << endl;
//     Array<double> input(resource.getData(), Serializable::Mode::JSON);
//     delta = deltaPadding(input);
//     delta_old = delta;
//     // unsigned int position = system_size/2;
//     // delta_start = delta[{position,position}];
// }

void Calculation::readDeltaCsv(int nr_sc_loop, string filename = "")
{
    if( filename == "")
    {
        filename = DeltaOutputFilename(nr_sc_loop, filename) + ".csv";
    }

    Array<double> input = DeltaFromCsv(filename);
    delta = deltaPadding(input);
    delta_old = delta;
    // unsigned int position = system_size/2;
    // delta_start = delta[{position,position}];
}


string Calculation::DeltaOutputFilename(const int &nr_sc_loop, const string& filename = "" )
{
    string nr_padded = to_string(nr_sc_loop);
    nr_padded.insert(0, 3 - nr_padded.length(), '0');
    if(filename != "")
    {
        return filename + "_delta_" + nr_padded;
    }
    else
    {
        return outputFileName + "_delta_" + nr_padded;
    }
    
}

Array<double> Calculation::realArray(const Array<complex<double>>& input)
{
    Array<double> out = Array<double>({input.getRanges()[0], input.getRanges()[1]}, 0);
    for(unsigned int i=0; i < input.getRanges()[0]; i++)
    {
        for(unsigned int j=0; j < input.getRanges()[1]; j++)
        {
            out[{i,j}] = real(input[{i,j}]);
        }
    }
    return out;
}

Array<complex<double>> Calculation::deltaPadding(const Array<double>& input)
{
    unsigned int sizeX_input = input.getRanges()[0];
    unsigned int sizeY_input = input.getRanges()[1];
    Array<complex<double>> out = Array<complex<double>>({system_size, system_size}, delta_start);
    if(system_size >= delta_simulation_size){
        unsigned int marginX = (system_size - sizeX_input)/2;
        unsigned int marginY = (system_size - sizeY_input)/2;
        for(unsigned int i=0; i < sizeX_input; i++)
        {
            for(unsigned int j=0; j < sizeY_input; j++)
            {
                out[{i+marginX,j+marginY}] = input[j+i*sizeY_input];
            }
        }
    }
    else{
        unsigned int marginX = (sizeX_input - system_size)/2;
        unsigned int marginY = (sizeY_input - system_size)/2;
        for(unsigned int i=0; i < system_size; i++)
        {
            for(unsigned int j=0; j < system_size; j++)
            {
                out[{i,j}] = input[(j+marginY)+(i+marginX)*sizeY_input];
            }
        }
    }

    return out;
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
    model = Model();
    impurity_level_present = false; //Will be set to true in functions that add an extra impurity layer
    for(int x = 0; x < system_size; x++){
        for(int y = 0; y < system_size; y++){
            for(unsigned int spin = 0; spin < 2; spin++){
                for(unsigned int ph = 0; ph < 2; ph++){
                    model << HoppingAmplitude(
                        -mu*(1. - 2*ph),
                        {x, y, spin, ph, 0},
                        {x, y, spin, ph, 0}
                    );
                    if(x+1 < system_size){
                        model << HoppingAmplitude(
                            -t*(1. - 2*ph),
                            {x+1, y, spin, ph, 0},
                            {x, y, spin, ph, 0}
                        ) + HC;
                    }
                    else
                    {
                        model << HoppingAmplitude(
                            -t*(1. - 2*ph),
                            {0, y, spin, ph, 0},
                            {x, y, spin, ph, 0}
                        ) + HC;
                    }
                    if(y+1 < system_size){
                        model << HoppingAmplitude(
                            -t*(1. - 2*ph),
                            {x, y+1, spin, ph, 0},
                            {x, y, spin, ph, 0}
                        ) + HC;
                    }
                    else
                    {
                        model << HoppingAmplitude(
                            -t*(1. - 2*ph),
                            {x, 0, spin, ph, 0},
                            {x, y, spin, ph, 0}
                        ) + HC; 
                    }
                }
                model << HoppingAmplitude(
                    // delta[{x,y}]*(1. - 2*spin),
                    functionDelta,
                    {x, y, spin, 0, 0},
                    {x, y, (spin+1)%2, 1, 0}
                ) + HC;
            }
            
            
        }
    }
    unsigned int position = system_size/2;
    // addImpurityLevel({position,position});
    // addPWaveBond({position,position});
    // addSOCBond({position,position});
    // addSwaveBond({position,position});
    // addSOC({position,position});
    // addLocalSwaveBonds({position,position});
    // addPWaveUP({position,position});
    // addPWave({position,position});
    for(unsigned int spin = 0; spin < 2; spin++){
        for(unsigned int ph = 0; ph < 2; ph++){
            model << HoppingAmplitude(
                Vz*(1. - 2*spin)*(1. - 2*ph),
                {position, position, spin, ph, 0},
                {position, position, spin, ph, 0}
            );
        }
    }
//     //Create model and set up hopping parameters
//     model = Model();

//     // Indeces: sample:0, tip: 1, x, y, spin
//     unsigned int system_index_sub = 0; //0 for the substrae, 1 for the tip
//     unsigned int system_index_imp;
//     unsigned int system_index_tip;

//     if(model_hubbard_model)
//     {
//         system_index_imp = system_index_sub + 1;
//         system_index_tip = system_index_sub + 2;
//     }
//     else
//     {       
//         system_index_imp = system_index_sub;
//         system_index_tip = system_index_sub + 1;
//     }
//     unsigned int position = system_size/2;
//     for(unsigned int x = 0; x < system_size; x++){
//         for(unsigned int y = 0; y < system_size; y++){
//             for(unsigned int s = 0; s < 2; s++){ // s... spin 0: up 1: down
// //-------------------BCS interaction term------------------------------------------

//                 model << HoppingAmplitude(Calculation::functionDelta, {system_index_sub, x, y, s, 0}, {system_index_sub, x, y, (s+1)%2, 1}) + HC;

//                 for(unsigned ph = 0; ph < 2; ph++){ // ph - particle hole, 0: particle, 1: hole
// // //------------------------chemical Potential-----------------------------------
//                     model << HoppingAmplitude(
//                         -mu*(1. - 2*ph),
//                         {system_index_sub, x, y, s, ph},
//                         {system_index_sub, x, y, s, ph}
//                     );
// //------------------------Nearest neighbour hopping term--------------------------------------
//                 //Add hopping parameters corresponding to t
//                     model << HoppingAmplitude(-t*(1. - 2*ph),	{system_index_sub,(x+1)%system_size, y, s, ph},	{system_index_sub,x, y, s, ph}) + HC;                
//                     model << HoppingAmplitude(-t*(1. - 2*ph),	{system_index_sub, x, (y+1)%system_size, s, ph},	{system_index_sub,x, y, s, ph}) + HC;

// //---------------------------Zeeman term------------------------------------------
//                 if(x == position and  y == position)
//                 {
//                     cout << Vz*(1. - 2*s)*(1. - 2*ph) << " @ " << system_index_imp << endl;
//                     model << HoppingAmplitude(Vz*(1. - 2*s)*(1. - 2*ph), {system_index_imp,x, y, s, ph}, {system_index_imp, x, y, s, ph});
//                 }
// //------------------------Rashba hopping term--------------------------------------
//                     model << HoppingAmplitude(-alpha*(1. - 2*ph)/2.,	{system_index_sub,x, y, 1, ph},	{system_index_sub,(x+1)%system_size, y, 0, ph}) + HC;
//                     model << HoppingAmplitude(alpha *(1. - 2*ph)/2.0,	{system_index_sub,(x+1)%system_size, y, 1, ph},	{system_index_sub,x, y, 0, ph}) + HC;
                
                
//                     model << HoppingAmplitude(-I*alpha*(1. - 2*ph)/2.,	{system_index_sub, x,y, 1, ph},	{system_index_sub,x,  (y+1)%system_size, 0, ph}) + HC;
//                     model << HoppingAmplitude(I*alpha*(1. - 2*ph)/2.,	{system_index_sub,x, (y+1)%system_size, 1, ph},	{system_index_sub,x, y, 0, ph}) + HC;
                  
//                 }
//             }
//         }
//     }





// //------------------------chemical Potential-----------------------------------
//                 //Add hopping amplitudes corresponding to chemical potential
//                 model << HoppingAmplitude(-mu,	{system_index_sub,x, y, s},	{system_index_sub,x, y, s});
//                 model << HoppingAmplitude(mu,	{system_index_sub,x, y, s+2},	{system_index_sub,x, y, s+2});

// //-------------------BCS interaction term------------------------------------------

// //                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
//                 model << HoppingAmplitude(Calculation::functionDelta, {system_index_sub,x,y,s}, {system_index_sub,x,y,(3-s)}) + HC;

// //-------------------BCS interaction term impurity------------------------------------------

// //                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
//                 if(x == system_size/2 and  y == system_size/2 and model_hubbard_model)
//                 {
//                     model << HoppingAmplitude(Calculation::functionDeltaDelta, {system_index_sub,x,y,s}, {system_index_sub,x,y,(3-s)}) + HC;
//                 }
// //------------------------Nearest neighbour hopping term--------------------------------------
//                 //Add hopping parameters corresponding to t
//                 if(x == system_size - 1){
//                     model << HoppingAmplitude(-t,	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(t,	{system_index_sub,x, y, s+2},{system_index_sub,(x+1)%system_size, y, s+2}) + HC;
//                 }
//                 else
//                 {
//                     model << HoppingAmplitude(-t,	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(t,	{system_index_sub,(x+1)%system_size, y, s+2},{system_index_sub,x, y, s+2}) + HC;
//                 }
                
//                 if(y == system_size - 1){
//                     model << HoppingAmplitude(-t,	{system_index_sub, x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(t,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
//                 }
//                 else
//                 {
//                     model << HoppingAmplitude(-t,	{system_index_sub,x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(t,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
//                 }
// //------------------------Rashba hopping term--------------------------------------
//                 if(x == system_size - 1){
//                     model << HoppingAmplitude(alpha *2.0*(0.5-s),	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(-alpha *2.0*(0.5-s),	{system_index_sub,x, y, s+2},{system_index_sub,(x+1)%system_size, y, s+2}) + HC;
//                 }
//                 else
//                 {
//                     model << HoppingAmplitude(alpha *2.0*(0.5-s),	{system_index_sub,(x+1)%system_size, y, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(-alpha *2.0*(0.5-s),	{system_index_sub,(x+1)%system_size, y, s+2},{system_index_sub,x, y, s+2}) + HC;
//                 }
                
//                 if(y == system_size - 1){
//                     model << HoppingAmplitude(i*alpha,	{system_index_sub, x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
//                 }
//                 else
//                 {
//                     model << HoppingAmplitude(-t,	{system_index_sub,x, (y+1)%system_size, s},	{system_index_sub,x, y, s}) + HC;
//                     model << HoppingAmplitude(t,  {system_index_sub,x, y, s+2}, {system_index_sub,x, (y+1)%system_size, s+2}) + HC;
//                 }  
//                 if(x == system_size - 1){
//                     model << HoppingAmplitude(alpha *2.0*(0.5-s), {system_index_sub, (x+1)%SIZE_X,y,(s+1)%2}, {system_index_sub, x,y,s}) + HC;
//                     model << HoppingAmplitude(-alpha *2.0*(0.5-s), {system_index_sub, x,y,s+2}, {system_index_sub, (x+1)%SIZE_X,y,(s+1)%2+2}) + HC;
//                 }

//                 if(periodicBoundCond || y+1 < SIZE_Y){
//     //                    model.addHAAndHC(HoppingAmplitude(i*alpha*2.0*(0.5-s),	{x, (y+1)%SIZE_Y, s*2},	{x, y, s*2+1}));
//     //                    model.addHAAndHC(HoppingAmplitude(-i*alpha*2.0*(0.5-s),  {x, y, s*2}, {x, (y+1)%SIZE_Y, s*2+1}));
//                     model.addHAAndHC(HoppingAmplitude(i*alpha, {x,(y+1)%SIZE_Y,(s+1)%2}, {x,y,s}));
//                     model.addHAAndHC(HoppingAmplitude(-i*alpha, {x,y,s+2}, {x,(y+1)%SIZE_Y,(s+1)%2+2}));

//                 }
                

// //---------------------------Zeeman term------------------------------------------
//                 if(x == position and  y == position)
//                 {
//                     model << HoppingAmplitude(Vz*2.0*(0.5-s), {system_index_imp,x, y, s}, {system_index_imp, x, y, s});
//                     model << HoppingAmplitude(-Vz*2.0*(0.5-s), {system_index_imp,x, y, s+2}, {system_index_imp,x, y, s+2});
//                     if(model_hubbard_model)
//                     {
//                         model << HoppingAmplitude(-mu,	{system_index_imp,x, y, s},	{system_index_imp,x, y, s});
//                         model << HoppingAmplitude(mu,	{system_index_imp,x, y, s+2},	{system_index_imp,x, y, s+2});
//                         model << HoppingAmplitude(-t_sample_imp,	{system_index_sub, x, y, s},	{system_index_imp, x, y, s}) + HC;
//                         model << HoppingAmplitude(t_sample_imp,  {system_index_sub, x, y, s+2}, {system_index_imp, x, y, s+2}) + HC;
//                     }
//                 }
//             }
//         }
//     }
//     if(model_hubbard_model)
//     {
//         //TODO Hubbard potential and hopping to impurity is still missing here
//     }
//     if(model_tip)
//     {
//         unsigned int position = tip_position;
//         if(!flat_tip)
//         {
//             for(unsigned int s = 0; s < 2; s++){
//                 model << HoppingAmplitude(-t_probe_sample,	{system_index_sub, position, position, s},	{system_index_tip, position, position, s}) + HC;
//                 model << HoppingAmplitude(t_probe_sample,  {system_index_sub, position, position, s+2}, {system_index_tip, position, position, s+2}) + HC;
            
//                 for(unsigned pos = 1; pos < probe_length; pos++){
//                     if(pos+1 < probe_length){
//                         model << HoppingAmplitude(-t_probe,	{pos, position, position, s},	{pos+1, position, position, s}) + HC;
//                         model << HoppingAmplitude(t_probe,  {pos, position, position, s+2}, {pos+1, position, position, s+2}) + HC;
//                     }
//                     model << HoppingAmplitude(-mu,	{pos, position, position, s},	{pos, position, position, s});
//                     model << HoppingAmplitude(mu,	{pos, position, position, s+2},	{pos, position, position, s+2});

//                     // model << HoppingAmplitude(Calculation::functionDeltaProbe, {system_index_tip, pos,s}, {system_index_tip, pos,(3-s)}) + HC;
//                     model << HoppingAmplitude(Calculation::functionDeltaProbe, {pos, position, position,s}, {pos, position, position,(3-s)}) + HC;
//                 }
//             }
//         }
//         else
//         {    
//             for(unsigned int s = 0; s < 2; s++){
//                 model << HoppingAmplitude(-t_probe_sample,	{system_index_sub, position, position, s},	{system_index_tip, position, position, s}) + HC;
//                 model << HoppingAmplitude(t_probe_sample,  {system_index_sub, position, position, s+2}, {system_index_tip, position, position, s+2}) + HC;
//             }
//             for(unsigned int x = 0; x < system_size; x++){
//                 for(unsigned int y = 0; y < system_size; y++){
//                     for(unsigned int s = 0; s < 2; s++){

//         //------------------------chemical Potential-----------------------------------
//                         //Add hopping amplitudes corresponding to chemical potential
//                         model << HoppingAmplitude(-mu,	{system_index_tip,x, y, s},	{system_index_tip,x, y, s});
//                         model << HoppingAmplitude(mu,	{system_index_tip,x, y, s+2},	{system_index_tip,x, y, s+2});

//         //-------------------BCS interaction term------------------------------------------

//         //                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
//                         model << HoppingAmplitude(Calculation::functionDeltaProbe, {system_index_tip,x,y,s}, {system_index_tip,x,y,(3-s)}) + HC;


//         //------------------------Nearest neighbour hopping term--------------------------------------
//                         //Add hopping parameters corresponding to t
//                         if(x == system_size - 1){
//                             model << HoppingAmplitude(-t_probe,	{system_index_tip,(x+1)%system_size, y, s},	{system_index_tip,x, y, s}) + HC;
//                             model << HoppingAmplitude(t_probe,	{system_index_tip,x, y, s+2},{system_index_tip,(x+1)%system_size, y, s+2}) + HC;
//                         }
//                         else
//                         {
//                             model << HoppingAmplitude(-t_probe,	{system_index_tip,(x+1)%system_size, y, s},	{system_index_tip,x, y, s}) + HC;
//                             model << HoppingAmplitude(t_probe,	{system_index_tip,(x+1)%system_size, y, s+2},{system_index_tip,x, y, s+2}) + HC;
//                         }
                        
//                         if(y == system_size - 1){
//                             model << HoppingAmplitude(-t_probe,	{system_index_tip, x, (y+1)%system_size, s},	{system_index_tip,x, y, s}) + HC;
//                             model << HoppingAmplitude(t_probe,  {system_index_tip,x, y, s+2}, {system_index_tip,x, (y+1)%system_size, s+2}) + HC;
//                         }
//                         else
//                         {
//                             model << HoppingAmplitude(-t_probe,	{system_index_tip,x, (y+1)%system_size, s},	{system_index_tip,x, y, s}) + HC;
//                             model << HoppingAmplitude(t_probe,  {system_index_tip,x, y, s+2}, {system_index_tip,x, (y+1)%system_size, s+2}) + HC;
//                         }

//                     }
//                 }
//             }
//         }
//     }
}

void Calculation::addImpurityLevel(const Index& pos){ //Impurity index is 1 (last position in Index)
    impurity_level_present = true;
    int x = pos.at(0);
    int y = pos.at(1);

        for(unsigned int spin = 0; spin < 2; spin++){
        for(unsigned int ph = 0; ph < 2; ph++){
            model << HoppingAmplitude(
                -mu*(1. - 2*ph),
                {x, y, spin, ph, 1},
                {x, y, spin, ph, 1}
            );
            model << HoppingAmplitude(
                -t_sample_imp*(1. - 2*ph),
                {x, y, spin, ph, 0},
                {x, y, spin, ph, 1}
            ) + HC;
            model << HoppingAmplitude(
                Vz*(1. - 2*spin)*(1. - 2*ph),
                {x, y, spin, ph, 1},
                {x, y, spin, ph, 1}
            );
        }
    }
}

void Calculation::addSOC(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    //Particle part x
    model << HoppingAmplitude(-alpha/2.,	{(x+1)%system_size, y, 1, 0, 0},	{x, y, 0, 0, 0});
    if(x-1 >= 0){
        model << HoppingAmplitude(alpha/2.,	{x-1, y, 1, 0, 0},	{x, y, 0, 0, 0});
        model << HoppingAmplitude(alpha/2.,	{x, y, 0, 0, 0},	{x-1, y, 1, 0, 0});
    }
    else{
        model << HoppingAmplitude(alpha/2.,	{system_size-1, y, 1, 0, 0},	{0, y, 0, 0, 0});
        model << HoppingAmplitude(alpha/2.,	{x, y, 0, 0, 0},	{system_size-1, y, 1, 0, 0});
    }
    model << HoppingAmplitude(-alpha/2.,	{x, y, 0, 0, 0},	{(x+1)%system_size, y, 1, 0, 0});

    //Particle part y
    model << HoppingAmplitude(-I*alpha/2.,	{x, (y+1)%system_size, 1, 0, 0},	{x, y, 0, 0, 0});
    if(y-1 >= 0){
        model << HoppingAmplitude(I*alpha/2.,	{x, y-1, 1, 0, 0},	{x, y, 0, 0, 0});
        model << HoppingAmplitude(-I*alpha/2.,	{x, y, 0, 0, 0},	{x, y-1, 1, 0, 0});
    }
    else{
        model << HoppingAmplitude(I*alpha/2.,	{x, system_size-1, 1, 0, 0},	{x, 0, 0, 0, 0});
        model << HoppingAmplitude(-I*alpha/2.,	{x, 0, 0, 0, 0},	{x, system_size-1, 1, 0, 0});
    }
    model << HoppingAmplitude(+I*alpha/2.,	{x, y, 0, 0, 0},	{x, (y+1)%system_size, 1, 0, 0});

    //Hole part x
    model << HoppingAmplitude(alpha/2.,	{x, y, 0, 1, 0},	{(x+1)%system_size, y, 1, 1, 0});
    if(x-1 >= 0){
        model << HoppingAmplitude(-alpha/2.,	{x, y, 0, 1, 0},	{(x-1), y, 1, 1, 0});
        model << HoppingAmplitude(-alpha/2.,	{x-1, y, 1, 1, 0},	{x, y, 0, 1, 0});
    }
    else{
        model << HoppingAmplitude(-alpha/2.,	{0, y, 0, 1, 0},	{system_size-1, y, 1, 1, 0});
        model << HoppingAmplitude(-alpha/2.,	{system_size-1, y, 1, 1, 0},	{0, y, 0, 1, 0});
    }
    model << HoppingAmplitude(alpha/2.,	{(x+1)%system_size, y, 1, 1, 0},	{x, y, 0, 1, 0});

    //Hole part y
    model << HoppingAmplitude(I*alpha/2.,	{x, y, 0, 1, 0},	{x, (y+1)%system_size, 1, 1, 0});
    if(y-1 >= 0){
        model << HoppingAmplitude(-I*alpha/2.,	{x, y, 0, 1, 0},	{x, y-1, 1, 1, 0});
        model << HoppingAmplitude(+I*alpha/2.,	{x, y-1, 1, 1, 0},	{x, y, 0, 1, 0});
    }
    else{
        model << HoppingAmplitude(-I*alpha/2.,	{x, 0, 0, 1, 0},	{x, system_size-1, 1, 1, 0});
        model << HoppingAmplitude(+I*alpha/2.,	{x, system_size-1, 1, 1, 0},	{x, 0, 0, 1, 0});
    }
    model << HoppingAmplitude(-I*alpha/2.,	{x, (y+1)%system_size, 1, 1, 0},	{x, y, 0, 1, 0});
}


void Calculation::addSOCBond(const Index& pos){
    impurity_level_present = true;
    int x = pos.at(0);
    int y = pos.at(1);

    //Particle part x
    model << HoppingAmplitude(-alpha/2.,	{x, y, 1, 0, 1},	{x, y, 0, 0, 0});
    model << HoppingAmplitude(-alpha/2.,	{x, y, 0, 0, 0},	{x, y, 1, 0, 1});
    model << HoppingAmplitude(alpha/2.,	{x, y, 1, 0, 0},	{x, y, 0, 0, 1});
    model << HoppingAmplitude(alpha/2.,	{x, y, 0, 0, 1},	{x, y, 1, 0, 0});


    //Hole part x
    model << HoppingAmplitude(alpha/2.,	{x, y, 0, 1, 0},	{x, y, 1, 1, 1});
    model << HoppingAmplitude(alpha/2.,	{x, y, 1, 1, 1},	{x, y, 0, 1, 0});
    model << HoppingAmplitude(-alpha/2.,	{x, y, 0, 1, 1},	{x, y, 1, 1, 0});
    model << HoppingAmplitude(-alpha/2.,	{x, y, 1, 1, 0},	{x, y, 0, 1, 1});
}


void Calculation::addPWaveBond(const Index& pos){
    impurity_level_present = true;
    int x = pos.at(0);
    int y = pos.at(1);
    for(unsigned spin = 0; spin < 2; spin++){
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 0, 0},
            {x, y, spin, 1, 1}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 1, 0},
            {x, y, spin, 0, 1}
        ) + HC;
    }
}

void Calculation::addSwaveBond(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    model << HoppingAmplitude(
        delta_s_bond,
        {x, y, 0, 0, 0},
        {x, y, 1, 1, 1}
    ) + HC;
    model << HoppingAmplitude(
        -delta_s_bond,
        {x, y, 1, 0, 0},
        {x, y, 0, 1, 1}
    ) + HC;
    model << HoppingAmplitude(
        -delta_s_bond,
        {x, y, 0, 1, 0},
        {x, y, 1, 0, 1}
    ) + HC;
    model << HoppingAmplitude(
        delta_s_bond,
        {x, y, 1, 1, 0},
        {x, y, 0, 0, 1}
    ) + HC;
}

void Calculation::addLocalSwaveBonds(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    for(int bond = 1; bond <= 1; bond += 2){
        // x direction
        model << HoppingAmplitude(
            delta_s_bond,
            {x, y, 0, 0, 0},
            {x + bond, y, 1, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            -delta_s_bond,
            {x, y, 1, 0, 0},
            {x + bond, y, 0, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            -delta_s_bond,
            {x, y, 0, 1, 0},
            {x+bond, y, 1, 0, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_s_bond,
            {x, y, 1, 1, 0},
            {x+bond, y, 0, 0, 0}
        ) + HC;
        // y direction
            model << HoppingAmplitude(
            delta_s_bond,
            {x, y, 0, 0, 0},
            {x, y + bond, 1, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            -delta_s_bond,
            {x, y, 1, 0, 0},
            {x, y + bond, 0, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            -delta_s_bond,
            {x, y, 0, 1, 0},
            {x, y+ bond, 1, 0, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_s_bond,
            {x, y, 1, 1, 0},
            {x, y + bond, 0, 0, 0}
        ) + HC;
    }

}

void Calculation::addUpDownPwaveBond(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    model << HoppingAmplitude(
        delta_p,
        {x, y, 0, 0, 0},
        {x, y, 1, 1, 1}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, 1, 0, 0},
        {x, y, 0, 1, 1}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, 0, 1, 0},
        {x, y, 1, 0, 1}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, 1, 1, 0},
        {x, y, 0, 0, 1}
    ) + HC;
}

void Calculation::addPWave(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    for(unsigned spin = 0; spin < 2; spin++){
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 0, 0},
            {x+1, y, spin, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 1, 0},
            {x+1, y, spin, 0, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 0, 0},
            {x-1, y, spin, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 1, 0},
            {x-1, y, spin, 0, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 0, 0},
            {x, y+1, spin, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 1, 0},
            {x, y+1, spin, 0, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 0, 0},
            {x, y-1, spin, 1, 0}
        ) + HC;
        model << HoppingAmplitude(
            delta_p,
            {x, y, spin, 1, 0},
            {x, y-1, spin, 0, 0}
        ) + HC;
    }

}

void Calculation::addPWaveUP(const Index& pos){
    int x = pos.at(0);
    int y = pos.at(1);

    unsigned spin = 0;

    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 0, 0},
        {x+1, y, spin, 1, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 1, 0},
        {x+1, y, spin, 0, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 0, 0},
        {x-1, y, spin, 1, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 1, 0},
        {x-1, y, spin, 0, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 0, 0},
        {x, y+1, spin, 1, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 1, 0},
        {x, y+1, spin, 0, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 0, 0},
        {x, y-1, spin, 1, 0}
    ) + HC;
    model << HoppingAmplitude(
        delta_p,
        {x, y, spin, 1, 0},
        {x, y-1, spin, 0, 0}
    ) + HC;

}


complex<double> Calculation::FunctionDelta::getHoppingAmplitude(const Index& from, const Index& to) const
{
    unsigned int from_x = from.at(0);
    unsigned int from_y = from.at(1);
    unsigned int from_s = from.at(2);
    unsigned int from_ph = from.at(3);
    // unsigned int from_ph = from.at(4);
    if(from_ph){
        return conj(delta[{from_x, from_y}])*(1. - 2*from_s)*(1. - 2*from_ph);
    } 
    else{
        return delta[{from_x, from_y}]*(1. - 2*from_s)*(1. - 2*from_ph);
    }

    // switch(from_s)
    // {
    // case 0:
    //     return conj(delta[{from_x, from_y}]);
    // case 1:
    //     return -conj(delta[{from_x, from_y}]);
    // case 2:
    //     return -delta[{from_x, from_y}];
    // case 3:
    //     return delta[{from_x, from_y}];
    // default:
    //     Streams::err << "something went wrong in Calculation::FuncDelta." << endl;
    //     return 0;
    // }
}

// complex<double> Calculation::FunctionDeltaDelta::getHoppingAmplitude(const Index& from, const Index& to) const
// {
//     unsigned int from_s = from.at(3);
//     return delta_probe*(1. - 2*from_s);
//     // switch(from_s)
//     // {
//     // case 0:
//     //     return -conj(delta_probe);
//     // case 1:
//     //     return conj(delta_probe);
//     // case 2:
//     //     return delta_probe;
//     // case 3:
//     //     return -delta_probe;
//     // default:
//     //     Streams::err << "something went wrong in Calculation::FuncDelta." << endl;
//     //     return 0;
//     // }
// }

// complex<double> Calculation::FunctionDeltaProbe::getHoppingAmplitude(const Index& from, const Index& to) const
// {
//     unsigned int from_s = from.at(3);
//     switch(from_s)
//     {
//     case 0:
//         return conj(delta_probe);
//     case 1:
//         return -conj(delta_probe);
//     case 2:
//         return -delta_probe;
//     case 3:
//         return delta_probe;
//     default:
//         Streams::err << "something went wrong in Calculation::FuncDelta." << endl;
//         return 0;
//     }
// }

#ifdef GPU_CALCULATION
bool Calculation::SelfConsistencyCallback::selfConsistencyCallback(Solver::ChebyshevExpander &solver)
#else
bool Calculation::SelfConsistencyCallback::selfConsistencyCallback(Solver::Diagonalizer &solver)
#endif 
{
    #ifdef GPU_CALCULATION
    PropertyExtractor::ChebyshevExpander pe(solver);
    #else
    PropertyExtractor::Diagonalizer pe(solver);
    #endif 

    pe.setEnergyWindow(-1*energy_bandwidth, 0, energy_points/2);






    delta_old = delta;
    Array<complex<double>> delta_temp = delta;
    double diff = 0.0;

    unsigned int position = system_size/2;
    for(unsigned int x=position; x < system_size; x++)
    {

        #pragma omp parallel for
        for(unsigned int y = position; y <= x; y++)
        {                   
            delta_temp[{x , y}] = (-pe.calculateExpectationValue({x,y, 0,1, 0},{x, y, 1,0,0})*coupling_potential*0.5 + delta_old[{x , y}]*0.5);
            if(abs((delta_temp[{x , y}]-delta_old[{x , y}]))/abs(delta_start) > diff)
            {
                diff = abs(delta_temp[{x , y}]-delta_old[{x , y}]);
            }
        }
    }
    diff = diff/abs(delta_start);

    for(unsigned int x=position; x < system_size; x++)
    {
        for(unsigned int y = position; y <= x; y++)
        {
            //Upper half
            //right
            delta[{x , y}] = delta_temp[{x , y}];
            delta[{y , x}] = delta_temp[{x , y}];
            //left
            delta[{2*position-x , y}] = delta_temp[{x , y}];
            delta[{y , 2*position-x}] = delta_temp[{x , y}];
            //Lower half
            //left
            delta[{2*position-x , 2*position-y}] = delta_temp[{x , y}];
            delta[{2*position-y , 2*position-x}] = delta_temp[{x , y}];
            //right
            delta[{x , 2*position-y}] = delta_temp[{x , y}];
            delta[{2*position-y , x}] = delta_temp[{x , y}];
        }
    }


    Streams::out << "Updated delta = " << to_string(real(delta_temp[{position,position}]/delta_start)) << ", ddelta = " << to_string(real(diff/delta_start)) << endl;
    if(abs(diff/delta_start) < EPS)
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
    #ifdef GPU_CALCULATION
    solver.setScaleFactor(4*energy_bandwidth);
    solver.setNumCoefficients(chebychev_coefficients);
    solver.setUseLookupTable(use_gpu);
    solver.setCalculateCoefficientsOnGPU(use_gpu);
    #endif
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
    Asolver.setModel(model);
    Asolver.setNumLanczosVectors(num_lanczos_vecs);
    Asolver.setMaxIterations(max_arnoldi_iterations);
    Asolver.setNumEigenValues(num_eigenvals);
    Asolver.setCalculateEigenVectors(false);
    Asolver.setCentralValue(-0.005);
    Asolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
    Asolver.run();
	Streams::out << "finished calc" << endl;
}

#ifndef GPU_CALCULATION
void Calculation::DoTestCalc()
{
    model.construct();
    solver.setModel(model);
    solver.run();
	Streams::out << "finished calc" << endl;
      //Extract DOS and write to file
    PropertyExtractor::Diagonalizer pe(solver);
	Property::DOS dos = pe.calculateDOS();
    Property::EigenValues ev = pe.getEigenValues();
    Exporter exporter;
    exporter.save(dos, outputFileName + "_dos.csv" );
    exporter.save(ev, outputFileName + "_Eigenvalues.csv");


    // int nr_excited_states = 20;
    // int middle = system_size*system_size*2;
    // vector<string> spin = {"e_up", "e_down", "h_up", "h_down"};


    
    // for(int j = middle - nr_excited_states/2; j < middle + nr_excited_states/2; j++){
    //   for(unsigned s = 0; s < 2; s++)
    //   {
    //     for(unsigned ph = 0; ph<2; ph++){
    //         Property::WaveFunctions wf = pe.calculateWaveFunctions(
    //         //   {{0, IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
    //         {{IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
    //         //   {system_size*4+i-nr_excited_states/2}
    //         {j}
    //          );

        
    //         ofstream output_real;
    //         output_real.open(outputFileName + "WaveFunction_" + spin[s + 2*ph] + "_nr_real_" + to_string(j) + ".csv", ios::trunc);
    //         ofstream output_imag;
    //         output_imag.open(outputFileName + "WaveFunction_" + spin[s + 2*ph] + "_nr_imag_" + to_string(j) + ".csv", ios::trunc);
    //         for(unsigned y = 0; y < system_size; y++)
    //         {
    //             for(unsigned x = 0; x < system_size; x++)
    //             {
    //                 output_real << real(wf({x,y,s,ph}, j));
    //                 output_imag << imag(wf({x,y,s,ph}, j));
    //                 if(x < system_size - 1)
    //                 {
    //                     output_real << ",";
    //                     output_imag << ",";
    //                 }   
    //             }
    //             output_real << endl;
    //             output_imag << endl;
    //         }
    //         output_real.close();
    //         output_imag.close();
    //     }
    //   }
	// }
    const double LOWER_BOUND = -5;
    const double UPPER_BOUND = 5;
    const unsigned int RESOLUTION = 500;
    pe.setEnergyWindow(
        LOWER_BOUND,
        UPPER_BOUND,
        RESOLUTION
    );
    Property::LDOS ldos = pe.calculateLDOS({
        {_a_, _a_, IDX_SUM_ALL, IDX_SUM_ALL}
    });

    const double SMOOTHING_SIGMA = 0.025;
    const unsigned int SMOOTHING_WINDOW = 51;
    ldos = Smooth::gaussian(ldos, SMOOTHING_SIGMA, SMOOTHING_WINDOW);





    ofstream ldos_file;
    ldos_file.open(outputFileName + "ldos_axisX.csv", ios::trunc);
    for(unsigned x = 0; x < system_size; x++){
        for(unsigned n = 0; n < ldos.getResolution(); n++){
            ldos_file << ldos({x,system_size/2,IDX_SUM_ALL,IDX_SUM_ALL}, n);
            if(n != ldos.getResolution() - 1){
                ldos_file << ",";
            }
        }
            if(x != system_size - 1){
                ldos_file << endl;
            }
    }
    ldos_file.close();
    ldos_file.open(outputFileName + "ldos_axisY.csv", ios::trunc);
    for(unsigned x = 0; x < system_size; x++){
        for(unsigned n = 0; n < ldos.getResolution(); n++){
            ldos_file << ldos({system_size/2,x,IDX_SUM_ALL,IDX_SUM_ALL}, n);
            if(n != ldos.getResolution() - 1){
                ldos_file << ",";
            }
        }
            if(x != system_size - 1){
                ldos_file << endl;
            }
    }
    ldos_file.close();
    ldos_file.open(outputFileName + "ldos_axisDiag.csv", ios::trunc);
    for(unsigned x = 0; x < system_size; x++){
        for(unsigned n = 0; n < ldos.getResolution(); n++){
            ldos_file << ldos({x,x,IDX_SUM_ALL,IDX_SUM_ALL}, n);
            if(n != ldos.getResolution() - 1){
                ldos_file << ",";
            }
        }
            if(x != system_size - 1){
                ldos_file << endl;
            }
    }
    ldos_file.close();
    // Extract LDOS and write to file
}
#endif

void Calculation::CalcEigenstates()
{
    int nr_excited_states = 5;
    double central_value = 0.01;

    model.construct();
    Asolver.setModel(model);
    Asolver.setNumLanczosVectors(30*nr_excited_states);
    Asolver.setMaxIterations(max_arnoldi_iterations*10);
    Asolver.setNumEigenValues(nr_excited_states);
    Asolver.setCalculateEigenVectors(calculate_waveFcts);
    Asolver.setCentralValue(-central_value);
    Asolver.setMode(Solver::ArnoldiIterator::Mode::ShiftAndInvert);
    Asolver.run();
	

    PropertyExtractor::ArnoldiIterator pe(Asolver);

	// Extract eigen values and write these to file
	Property::EigenValues ev = pe.getEigenValues();
    Exporter exporter;
    exporter.save(ev, outputFileName + "Eigenvalues_minus.csv" );


    vector<string> spin = {"e_up", "e_down", "h_up", "h_down"};

    if(calculate_waveFcts){
        Property::WaveFunctions wf = pe.calculateWaveFunctions(
        //   {{0, IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
            {{IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
        //   {system_size*4+i-nr_excited_states/2}
            {IDX_ALL}
        );



        for(int j = 0; j < nr_excited_states; j++){
            for(unsigned s = 0; s < 2; s++)
            {
                for(unsigned ph = 0; ph<2; ph++){

                
                    ofstream output_real;
                    output_real.open(outputFileName + "WaveFunction_minus_" + spin[s + 2*ph] + "_nr_real_" + to_string(j) + ".csv", ios::trunc);
                    ofstream output_imag;
                    output_imag.open(outputFileName + "WaveFunction_minus_" + spin[s + 2*ph] + "_nr_imag_" + to_string(j) + ".csv", ios::trunc);
                    for(unsigned y = 0; y < system_size; y++)
                    {
                        for(unsigned x = 0; x < system_size; x++)
                        {
                            // output_real << real(wf({0,x,y,s,ph}, j));
                            // output_imag << imag(wf({0,x,y,s,ph}, j));
                            output_real << real(wf({x,y,s,ph, 0}, j));
                            output_imag << imag(wf({x,y,s,ph, 0}, j));
                            if(x < system_size - 1)
                            {
                                output_real << ",";
                                output_imag << ",";
                            }   
                        }
                        output_real << endl;
                        output_imag << endl;
                    }
                    output_real.close();
                    output_imag.close();
                    if(impurity_level_present){
                        unsigned middle = system_size/2;
                        output_real.open(outputFileName + "WaveFunction_minus_imp_" + spin[s + 2*ph] + "_nr_real_" + to_string(j) + ".csv", ios::trunc);
                        output_imag.open(outputFileName + "WaveFunction_minus_imp" + spin[s + 2*ph] + "_nr_imag_" + to_string(j) + ".csv", ios::trunc);
                        output_real << real(wf({middle, middle,s,ph, 1}, j));
                        output_imag << imag(wf({middle, middle,s,ph, 1}, j));
                        output_real << endl;
                        output_imag << endl;
                        output_real.close();
                        output_imag.close();
                    }
                }
            }
        }
    }
    
    Asolver.setCentralValue(central_value);
    Asolver.run();
    ev = pe.getEigenValues();
    exporter.save(ev, outputFileName + "Eigenvalues_plus.csv" );


    if(calculate_waveFcts){
        Property::WaveFunctions wf = pe.calculateWaveFunctions(
            //   {{0, IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
            {{IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
            //   {system_size*4+i-nr_excited_states/2}
            {IDX_ALL}
        );

        for(int j = 0; j < nr_excited_states; j++){
            for(unsigned s = 0; s < 2; s++)
            {
                for(unsigned ph = 0; ph<2; ph++){

                
                    ofstream output_real;
                    output_real.open(outputFileName + "WaveFunction_plus_" + spin[s + 2*ph] + "_nr_real_" + to_string(j) + ".csv", ios::trunc);
                    ofstream output_imag;
                    output_imag.open(outputFileName + "WaveFunction_plus_" + spin[s + 2*ph] + "_nr_imag_" + to_string(j) + ".csv", ios::trunc);
                    for(unsigned y = 0; y < system_size; y++)
                    {
                        for(unsigned x = 0; x < system_size; x++)
                        {
                            // output_real << real(wf({0,x,y,s,ph}, j));
                            // output_imag << imag(wf({0,x,y,s,ph}, j));
                            output_real << real(wf({x,y,s,ph, 0}, j));
                            output_imag << imag(wf({x,y,s,ph, 0}, j));
                            if(x < system_size - 1)
                            {
                                output_real << ",";
                                output_imag << ",";
                            }   
                        }
                        output_real << endl;
                        output_imag << endl;
                    }
                    output_real.close();
                    output_imag.close();
                    if(impurity_level_present){
                        unsigned middle = system_size/2;
                        output_real.open(outputFileName + "WaveFunction_plus_imp_" + spin[s + 2*ph] + "_nr_real_" + to_string(j) + ".csv", ios::trunc);
                        output_imag.open(outputFileName + "WaveFunction_plus_imp" + spin[s + 2*ph] + "_nr_imag_" + to_string(j) + ".csv", ios::trunc);
                        output_real << real(wf({middle, middle,s,ph, 1}, j));
                        output_imag << imag(wf({middle, middle,s,ph, 1}, j));
                        output_real << endl;
                        output_imag << endl;
                        output_real.close();
                        output_imag.close();
                    }

                }
            }
        }
    }
    Streams::out << "finished calc" << endl;
}

void Calculation::WriteOutputSc()
{
    // PropertyExtractor::ChebyshevExpander pe(solver);
    // Exporter exporter;
    // exporter.setFormat(Exporter::Format::ColumnMajor);
    #ifndef GPU_CALCULATION 
    PropertyExtractor::Diagonalizer pe(solver);
    //Extract eigen values and write these to file
    Property::EigenValues ev = pe.getEigenValues();
    Exporter exporter;
    exporter.save(ev, outputFileName + "Eigenvalues.csv" );
    #endif
    // FileWriter::setFileName(outputFileName);

    // const double UPPER_BOUND = 5; //10*abs(delta_start);
	// const double LOWER_BOUND = -5; //-10*abs(delta_start);
	// const int RESOLUTION = 2000;
	// pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, RESOLUTION);


    // FileWriter::setFileName(outputFileName);

  //Extract DOS and write to file
	// Property::DOS dos = pe.calculateDOS();
	// FileWriter::writeDOS(dos);



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



//   int nr_excited_states = 30;

//   for(int i = 1; i <= nr_excited_states; i++){
//       Property::WaveFunctions wf = pe.calculateWaveFunctions(
//           {{IDX_ALL, IDX_ALL, IDX_ALL}},
//           {system_size*4+i-nr_excited_states/2}
//       );
//       FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));

//   }

}


void Calculation::WriteOutput()
{
    PropertyExtractor::ArnoldiIterator pe(Asolver);
    // PropertyExtractor::Diagonalizer pe(solver);
    // FileWriter::setFileName(outputFileName);

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

	// Extract eigen values and write these to file
	Property::EigenValues ev = pe.getEigenValues();
    Exporter exporter;
    exporter.save(ev, outputFileName + "Eigenvalues.csv" );

	// Extract LDOS and write to file

        // Property::LDOS ldos = pe.calculateLDOS(
        //     {IDX_Z, IDX_X, IDX_Y, IDX_SUM_ALL},
        //     {1, system_size, system_size,	4}
        // );
        // FileWriter::writeLDOS(ldos);

	// Extract LDOS and write to file

        // Property::LDOS ldos = pe.calculateLDOS({
        //     {0, system_size/2, system_size/2, IDX_SUM_ALL},
        //     {1, system_size/2, system_size/2, IDX_SUM_ALL},
        // });
        // FileWriter::writeLDOS(ldos);

//     int nr_excited_states = 30;

//   for(int i = 1; i <= nr_excited_states; i++){
//       Property::WaveFunctions wf = pe.calculateWaveFunctions(
//           {{0, IDX_ALL, IDX_ALL, IDX_ALL}},
//           {system_size*4+i-nr_excited_states/2}
//       );
//       FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));
// 	}

//   int nr_excited_states = 150;

//   for(int i = 0; i < nr_excited_states; i++){
//       Property::WaveFunctions wf = pe.calculateWaveFunctions(
//           {{IDX_ALL, IDX_ALL, IDX_ALL, IDX_ALL}},
//           {i}
//       );
//       FileWriter::writeWaveFunctions(wf, "WaveFunction_" + to_string(i));

//   }


   WriteDelta(0);
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
    string filename = DeltaOutputFilename(nr_loop);
    Exporter exporter;
    exporter.save(GetRealVec(delta), filename + ".csv" );

    // if(nr_loop == 0)
    // {
    //     Array<double> delta_real = realArray(delta);
    //     string delta_out = delta_real.serialize(Serializable::Mode::JSON);
    //     Resource resource;
    //     resource.setData(delta_out);
    //     resource.write(filename + ".json");
    // }

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

Array<complex<double>> Calculation::CalculateChargeDensity(unsigned int spin)
{
    #ifdef GPU_CALCULATION
	PropertyExtractor::ChebyshevExpander pe(solver);
    #else
	PropertyExtractor::Diagonalizer pe(solver);
    #endif
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

void Calculation::setAlpha(complex<double> input)
{
  alpha = input;
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

void Calculation::setTipPosition(unsigned int position)
{
    tip_position = position;
}

unsigned int Calculation::getSystemSize()
{
    return system_size;
}

complex<double> Calculation::getDeltaStart()
{
    return delta_start;
}

void Calculation::setDeltaDelta(complex<double> dD)
{
    delta_Delta = dD;
}

Array<double> Calculation::ConvertVectorToArray(const double *input, unsigned int sizeX, unsigned int sizeY)
{
    Array<double> out = Array<double>({sizeX, sizeX}, 0);
    
    for(unsigned int i=0; i < sizeX; i++)
    {
        for(unsigned int j=0; j < sizeY; j++)
        {
            out[{i,j}] = input[j+i*sizeY];
        }
    }
    return out;
}


Array<double> Calculation::DeltaFromCsv(const string& filename)
{
    vector<double> temp_vec;
    std::ifstream csv_file(filename);
    string line;

    while(csv_file.good())
    {
        getline(csv_file, line);
        try
        {
            temp_vec.push_back(stod(line));
        }
        catch(const std::invalid_argument& e)
        {
            std::cerr << line << '\n';
        }
    }
    if(temp_vec.size() != delta_simulation_size*delta_simulation_size)
    {
        TBTKExit("Calculation::DeltaFromCsv", "Wrong size of delta input file", "Check if datapoints agree with size " + to_string(delta_simulation_size));
    }

    return ConvertVectorToArray(temp_vec.data(), delta_simulation_size, delta_simulation_size);
}
