#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>

#include "TBTK/Model.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Solver/ChebyshevExpander.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/PropertyExtractor/ChebyshevExpander.h"
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
        void WriteOutput();
        void setVz(complex<double>);
        void setMu(complex<double>);
        void setOutputFileName(string);
        void setcoupling_potential(complex<double>);
        void AddDefects(int);



    protected:

    private:

        void WriteDelta(int);
        Array<double> GetRealVec(Array<complex<double>>);
        Array<double> GetImagVec(Array<complex<double>> );

        class FunctionDelta : 
        public HoppingAmplitude::AmplitudeCallback
        {
                complex<double> getHoppingAmplitude(const Index&, const Index&) const;    
            
        };
        static FunctionDelta functionDelta;
        
        // public Solver::Diagonalizer::SelfConsistencyCallback
        class SelfConsistencyCallback :
        public Solver::ChebyshevExpander::SelfConsistencyCallback
        {
                public:
                bool selfConsistencyCallback(Solver::ChebyshevExpander &solver);
                // bool selfConsistencyCallback(Solver::Diagonalizer &solver);
                
        };
        static SelfConsistencyCallback selfConsistencyCallback;

        static unsigned int system_length;
        static unsigned int system_size;
        static complex<double> mu;
        static complex<double> Vz;
        static complex<double> t;
        static complex<double> delta_start;
        static complex<double> coupling_potential;

        static const complex<double> I;
        static const double EPS;

        // static Solver::Diagonalizer solver;
        static Solver::ChebyshevExpander solver;
        static Array<complex<double>> delta;
        static Array<complex<double>> delta_old;

        static string outputFileName;

        bool symmetry_on;
        bool use_gpu;


        Model model;
};

#endif // CALCULATION_H
