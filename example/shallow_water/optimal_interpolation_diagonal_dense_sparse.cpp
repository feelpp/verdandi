// Program to test the equivalence between dense and sparse observation
// operator. Use the identity observation operator to do so either dense or
// sparse. Example to be called with river-sparse.cfg,
// river-diagonal-dense.cfg and river-operator-from-file.cfg.
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU
#define SELDON_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT
#define GETPOT_ACTIVATE_EXCEPTION true

#include "Verdandi.hxx"
using namespace Verdandi;

#include "seldon/SeldonSolver.hxx"

#include "OptimalInterpolation.cxx"
#include "LinearObservationManager.cxx"
#include "ShallowWater.cxx"
#include "newran.h"

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
    typedef ShallowWater<real> ClassModel;
    typedef OptimalInterpolation<real, ClassModel,
        LinearObservationManager<real> > ClassOptimalInterpolation;

    ClassOptimalInterpolation driver(argv[1]);

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();

        driver.Forward();

        driver.Analyze();
    }

    END;

    return 0;

}
