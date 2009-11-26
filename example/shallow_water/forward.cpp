// Generates the truth observations by applying the model from truth initial
// conditions (river-truth.cfg).
#define SELDON_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT
#define GETPOT_ACTIVATE_EXCEPTION true

#include "Verdandi.hxx"
using namespace Verdandi;

#include "ForwardDriver.cxx"

#include "ShallowWater.cxx"

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

    ForwardDriver<ShallowWater<double> > driver(argv[1]);

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    END;

    return 0;

}
