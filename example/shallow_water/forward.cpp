// Generates the truth observations by applying the model from truth initial
// conditions (river-truth.cfg).
#define SELDON_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT
#define GETPOT_ACTIVATE_EXCEPTION true

#include "Verdandi.hxx"
using namespace Verdandi;

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

    ShallowWater<double> shallow_water(argv[1]);
    shallow_water.Initialize(argv[1]);

    while (!shallow_water.HasFinished())
    {
        shallow_water.InitializeStep();

        shallow_water.Forward();
    }

    END;

    return 0;

}
