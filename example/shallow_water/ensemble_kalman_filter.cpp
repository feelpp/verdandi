#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK


#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"

#include "model/ShallowWater.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/EnsembleKalmanFilter.cxx"
#include "method/NewranPerturbationManager.cxx"

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

    typedef double real;

    Verdandi::EnsembleKalmanFilter<real, Verdandi::ShallowWater<real>,
        Verdandi::LinearObservationManager<real>,
        Verdandi::NewranPerturbationManager> driver;

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
