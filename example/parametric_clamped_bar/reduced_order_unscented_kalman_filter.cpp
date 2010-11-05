#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_WITH_ABORT
#define VERDANDI_SPARSE

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS
//#define VERDANDI_WITH_MPI
//#define VERDANDI_WITH_OMP

//#define VERDANDI_LOGGING_LEVEL -7

#if defined(VERDANDI_WITH_MPI)
#include <mpi.h>
#endif

#if defined(VERDANDI_WITH_OMP)
#include <omp.h>
#endif

#include "Verdandi.hxx"
using namespace Verdandi;
#include "seldon/SeldonSolver.hxx"

#include "model/ParametricClampedBar.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/ReducedOrderUnscentedKalmanFilter.cxx"


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

#if defined(VERDANDI_WITH_MPI)
    MPI::Init(argc, argv);
#endif

    typedef double real;

    ReducedOrderUnscentedKalmanFilter<real, ParametricClampedBar<real>,
        LinearObservationManager<real> > driver(argv[1]);

    driver.Initialize();

    while (!driver.HasFinished())
    {
        driver.InitializeStep();

        driver.Forward();

        driver.Analyze();
    }

#if defined(VERDANDI_WITH_MPI)
    MPI::Finalize();
#endif

    END;

    return 0;

}
