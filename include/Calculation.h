#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>

#include "TBTK/Model.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/Solver/ArnoldiIterator.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
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
        void setMu(complex<double>);
        void setPhase(double);
        void setOutputFileName(string);
        void setcoupling_potential(complex<double>);
        void AddDefects(int);
        void readDelta(int, string);
        void WriteDelta(int);



    protected:

    private:

        
        Array<double> GetRealVec(Array<complex<double>>);
        Array<double> GetImagVec(Array<complex<double>> );
        Array<complex<double>> CalculateChargeDensity(unsigned int);
        Array<complex<double>> ConvertVectorToArray(const double *, unsigned int, unsigned int);
        Array<complex<double>> deltaPadding(Array<complex<double>>, unsigned int, unsigned int, unsigned int, unsigned int);

        class FunctionDelta : 
        public HoppingAmplitude::AmplitudeCallback
        {
                complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        };
        static FunctionDelta functionDelta;

        class FunctionDeltaProbe : 
        public HoppingAmplitude::AmplitudeCallback
        {
                complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        };
        static FunctionDeltaProbe functionDeltaProbe;
        
        class SelfConsistencyCallback :
        public Solver::Diagonalizer::SelfConsistencyCallback
        {
                public:
                bool selfConsistencyCallback(Solver::Diagonalizer &solver);
                
        };
        static SelfConsistencyCallback selfConsistencyCallback;

        static unsigned int system_length;
        static unsigned int system_size;
        static complex<double> mu;
        static complex<double> mu_probe;
        static complex<double> Vz;
        static complex<double> t;
        static complex<double> t_probe;
        static complex<double> t_probe_sample;
        static complex<double> t_sample_imp;
        static double phase;
        static unsigned int probe_length;
        static complex<double> delta_start;
        static complex<double> delta_probe;
        static complex<double> coupling_potential;

        static const complex<double> I;
        static const double EPS;

        static Solver::Diagonalizer solver;
        static Solver::ArnoldiIterator Asolver;
        //static Solver::ChebyshevExpander solver;
        static Array<complex<double>> delta;
        static Array<complex<double>> delta_old;

        static string outputFileName;

        bool symmetry_on;
        bool use_gpu;
        bool model_tip;


        Model model;
};

#endif // CALCULATION_H
