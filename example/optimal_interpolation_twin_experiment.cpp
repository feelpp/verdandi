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

    typedef double real;
    typedef ShallowWater<real> ClassModel;
    typedef OptimalInterpolation<real, ClassModel,
	GridToNetworkObservationManager<real> > ClassOptimalInterpolation;

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
