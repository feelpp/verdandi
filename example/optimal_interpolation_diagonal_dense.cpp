// Program to test the equivalence between dense and sparse observation
// operator. Use the identity observation operator to do so either dense or
// sparse.
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_ABORT

#include "Verdandi.hxx"
using namespace Verdandi;

#include "OptimalInterpolation.cxx"
#include "LinearObservationManager.cxx"
#include "ShallowWater.cxx"
#include "OutputSaver.cxx"
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

    OutputSaver<real, ClassOptimalInterpolation> output_saver(argv[1],
                                                              driver);

    driver.Initialize(argv[1]);
    output_saver.Initialize(argv[1], driver);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        output_saver.InitializeStep();
        driver.Forward();
        output_saver.Save(driver);
        driver.Analyze();
        output_saver.Save(driver);
    }

    END;

    return 0;

}
