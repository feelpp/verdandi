#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"
using namespace Verdandi;

#include "model/ModelTemplate.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/OptimalInterpolation.cxx"


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

    OptimalInterpolation<double, ModelTemplate,
        LinearObservationManager<double> > driver(argv[1]);

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