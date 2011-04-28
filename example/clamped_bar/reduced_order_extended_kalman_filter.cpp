#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE

//#define VERDANDI_WITH_DIRECT_SOLVER
//#define SELDON_WITH_MUMPS

#if defined(VERDANDI_WITH_MPI)
#include <mpi.h>
#endif


#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"

#include "model/ClampedBar.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/ReducedOrderExtendedKalmanFilter.cxx"


int main(int argc, char** argv)
{

    TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        std::cout << mesg << std::endl;
        return 1;
    }

#if defined(VERDANDI_WITH_MPI)
    MPI::Init(argc, argv);
#endif

    typedef double real;

    Verdandi::ReducedOrderExtendedKalmanFilter<real,
        Verdandi::ClampedBar<real>,
        Verdandi::LinearObservationManager<real> > driver;

    driver.Initialize(argv[1]);

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
