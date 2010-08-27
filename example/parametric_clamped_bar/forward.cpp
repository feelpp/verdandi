#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE

#define VERDANDI_WITH_DIRECT_SOLVER
#define SELDON_WITH_SUPERLU

#include "Verdandi.hxx"
using namespace Verdandi;
#include "seldon/SeldonSolver.hxx"

#include "model/ParametricClampedBar.cxx"
#include "method/ForwardDriver.cxx"


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

    ForwardDriver<ParametricClampedBar<double> > driver(argv[1]);

    driver.Initialize();

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    END;

    return 0;

}
