// You may use 'LinearObservationManager' or 'GridToNetworkObservationManager'
// as observation operator.
#ifndef OBSERVATION_OPERATOR
#define OBSERVATION_OPERATOR LinearObservationManager
#endif


#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU
#define VERDANDI_SPARSE

#include "Verdandi.hxx"
using namespace Verdandi;

#include "seldon/SeldonSolver.hxx"

#include "OptimalInterpolation.cxx"
#include "ShallowWater.cxx"

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#include QUOTE(OBSERVATION_OPERATOR.cxx)


int main(int argc, char** argv)
{

    TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        cout << mesg << endl;
        return 1;
    }

    typedef double real;

    OptimalInterpolation<real, ShallowWater<real>,
        OBSERVATION_OPERATOR<real> > driver(argv[1]);

    driver.Initialize();

    while (!driver.HasFinished())
    {
        driver.InitializeStep();

        driver.Forward();

        driver.Analyze();
    }

    END;

    return 0;

}
