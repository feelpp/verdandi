#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK

//#define SELDON_WITH_MUMPS
#define SELDON_WITH_SUPERLU
//#define SELDON_WITH_UMFPACK

//#define SELDON_WITH_MPI

#define VERDANDI_WITH_DIRECT_SOLVER

#if !defined(SELDON_WITH_UMFPACK) && !defined(SELDON_WITH_SUPERLU)      \
    && !defined(SELDON_WITH_MUMPS) && defined(VERDANDI_WITH_DIRECT_SOLVER)
#define SELDON_WITH_SUPERLU
#endif

#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"

using namespace Verdandi;
#include "seldon/SeldonSolver.hxx"

#include "ExtendedKalmanFilter.cxx"
#include "LinearObservationManager.cxx"
#include "ClampedBar.cxx"


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

    ExtendedKalmanFilter<real, ClampedBar<real>,
        LinearObservationManager<real> > driver(argv[1]);

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
