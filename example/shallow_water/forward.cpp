#define SELDON_DEBUG_LEVEL_4

#define VERDANDI_WITH_ABORT

#include "Verdandi.hxx"

#include "model/ShallowWater.cxx"
#include "method/ForwardDriver.cxx"


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

    Verdandi::ForwardDriver<Verdandi::ShallowWater<double> > driver(argv[1]);

    driver.Initialize();

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    END;

    return 0;

}
