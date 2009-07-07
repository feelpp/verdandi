// Generates the truth observations by applying the model from truth initial
// conditions (river-truth.cfg).
#define SELDONDATA_DEBUG_LEVEL_4
#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_ABORT

#include <iostream>
using namespace std;

#include "OptimalInterpolation.cxx"
#include "GridToNetworkObservationManager.cxx"
#include "ShallowWater.cxx"
#include "OutputSaver.cxx"
#include "newran.h"

namespace SeldonData
{
    using Talos::to_num;
    using Talos::to_str;
    class Error;
    class IOError;
    class WrongDim;
}

using namespace Verdandi;

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

    OutputSaver<double, ShallowWater<double> > output_saver(argv[1],
                                                            shallow_water);
    output_saver.Initialize(argv[1], shallow_water);

    while (!shallow_water.HasFinished())
    {
        shallow_water.InitializeStep();
        output_saver.InitializeStep();
        shallow_water.Forward();
        output_saver.Save(shallow_water);
    }

    END;

    return 0;

}

