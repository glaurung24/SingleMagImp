#pragma once

#define GPU_CALCULATION

#include <complex>

#include "TBTK/Model.h"

#ifdef GPU_CALCULATION
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
#else
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#endif
#include "TBTK/Solver/ArnoldiIterator.h"
#include "TBTK/PropertyExtractor/ArnoldiIterator.h"
#include "TBTK/Array.h"



using namespace std;
using namespace TBTK;

class Calculation
{
    public:
        Calculation();
        Calculation(string, complex<double>);
        virtual ~Calculation();
        void Init(string outputfilename, complex<double> vz_input = 0.0);
        void InitModel();
        void DoScCalc();
        void DoCalc();
        void WriteOutput();
        void WriteOutputSc();
        void setVz(complex<double>);
        void setAlpha(complex<double>);
        void setTImpSample(complex<double>);
        void setMu(complex<double>);
        void setPhase(double);
        void setOutputFileName(string);
        void setcoupling_potential(complex<double>);
        void setTipPosition(unsigned int);
	void setSystem_length(unsigned int);
        void AddDefects(int);
        // void readDeltaHdf5(int, string);
        void readDeltaCsv(int, string);
        void readDeltaJson(int, string);
        void WriteDelta(int);
        unsigned int getSystemSize();
        complex<double> getDeltaStart();
        void setDeltaDelta(complex<double>);
        string DeltaOutputFilename(const int& nr_sc_loop, const string& filename);
	void runArnoldiIterator();
        void CalcEigenstates();
        void DoTestCalc();



    protected:

    private:

        
        Array<double> GetRealVec(Array<complex<double>>);
        Array<double> GetImagVec(Array<complex<double>>);
        Array<double> realArray(const Array<complex<double>>&);
        Array<double> DeltaFromCsv(const string&);
        Array<complex<double>> CalculateChargeDensity(unsigned int);
        Array<double> ConvertVectorToArray(const double *, unsigned int, unsigned int);
        Array<complex<double>> deltaPadding(Array<complex<double>>, unsigned int, unsigned int, unsigned int, unsigned int);
        Array<complex<double>> deltaPadding(const Array<double>&);
        void addSOC(const Index& pos);
        void addPWave(const Index& pos);
        void addImpurityLevel(const Index& pos);
        void addPWaveBond(const Index& pos);
        void addSOCBond(const Index& pos);
        void addSwaveBond(const Index& pos);
        void addUpDownPwaveBond(const Index& pos);
        void addLocalSwaveBonds(const Index& pos);
        void addPWaveUP(const Index& pos);
        

        class FunctionDelta : 
        public HoppingAmplitude::AmplitudeCallback
        {
                complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        };
        static FunctionDelta functionDelta;

        // class FunctionDeltaProbe : 
        // public HoppingAmplitude::AmplitudeCallback
        // {
        //         complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        // };
        // static FunctionDeltaProbe functionDeltaProbe;

        // class FunctionDeltaDelta : 
        // public HoppingAmplitude::AmplitudeCallback
        // {
        //         complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        // };
        // static FunctionDeltaDelta functionDeltaDelta;

	#ifdef GPU_CALCULATION
    class SelfConsistencyCallback
	{
		public:
                bool selfConsistencyCallback(Solver::ChebyshevExpander &solver);
        };
        static SelfConsistencyCallback selfConsistencyCallback;
	#else
    class SelfConsistencyCallback :
        public Solver::Diagonalizer::SelfConsistencyCallback
        {
		public:
                bool selfConsistencyCallback(Solver::Diagonalizer &solver);
                
        };
        static SelfConsistencyCallback selfConsistencyCallback;
	#endif
        

        static unsigned int system_length;
        static unsigned int delta_simulation_size;
        static unsigned int system_size;
        static unsigned int tip_position;
        static unsigned int energy_points;
        static unsigned int chebychev_coefficients;
        static double lower_energy_bound;
        static double upper_energy_bound;
        static unsigned int max_arnoldi_iterations;
        static unsigned int num_eigenvals;
        static unsigned int num_lanczos_vecs;
        static unsigned int max_sc_iterations;
        static double energy_bandwidth;
        static complex<double> mu;
        static complex<double> mu_probe;
        static complex<double> Vz;
        static complex<double> t;
        static complex<double> t_probe;
        static complex<double> t_probe_sample;
        static complex<double> t_sample_imp;
        static complex<double> alpha;
        static complex<double> delta_p;
        static complex<double> delta_s_bond;
        static double phase;
        static unsigned int probe_length;
        static complex<double> delta_start;
        static complex<double> delta_probe;
        static complex<double> delta_Delta;
        static complex<double> coupling_potential;

        static const complex<double> I;
        static const double EPS;

        static Solver::ArnoldiIterator Asolver;
        #ifdef GPU_CALCULATION
        static Solver::ChebyshevExpander solver;
        #else
        static Solver::Diagonalizer solver;
        #endif
        
        static Array<complex<double>> delta;
        static Array<complex<double>> delta_old;

        static string outputFileName;

        static bool symmetry_on;
        static bool use_gpu;
        static bool model_tip;
        static bool flat_tip;
        bool model_hubbard_model;
        bool calculate_waveFcts;
        bool impurity_level_present;


        Model model;
};
