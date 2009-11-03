// Generates the truth observations by applying the model from truth initial
// conditions (river-truth.cfg).
#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_ABORT
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK

//#define SELDON_WITH_MUMPS
//#define SELDON_WITH_SUPERLU
//#define SELDON_WITH_UMFPACK

//#define SELDON_WITH_MPI

#define VERDANDI_WITH_ABORT

#define VERDANDI_WITH_DIRECT_SOLVER

#if !defined(SELDON_WITH_UMFPACK) && !defined(SELDON_WITH_SUPERLU)      \
    && !defined(SELDON_WITH_MUMPS) && defined(VERDANDI_WITH_DIRECT_SOLVER)
//#define SELDON_WITH_UMFPACK
//#define SELDON_WITH_SUPERLU
#define SELDON_WITH_MUMPS
#endif



#include "Verdandi.hxx"
using namespace Verdandi;

#include "ClampedBar.cxx"
#include "newran/newran.h"

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

    ClampedBar<double> clamped_bar(argv[1]);
    clamped_bar.Initialize(argv[1]);
    clamped_bar.InitializeFirstStep();
    /*
      OutputSaver<double, ClampedBar<double> > output_saver(argv[1],
      clamped_bar);
      output_saver.Initialize(argv[1], clamped_bar);
    */
    while (!clamped_bar.HasFinished())
    {
        clamped_bar.InitializeStep();
        //output_saver.InitializeStep();
        clamped_bar.Forward();
        //output_saver.Save(clamped_bar);
    }

    END;

    return 0;

}


