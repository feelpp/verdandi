// Program to write a matrix in a file.
#define SELDON_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT
#define VERDANDI_IGNORE_MESSAGE

#include "Verdandi.hxx"
using namespace Verdandi;

int main(int argc, char** argv)
{
    TRY;

    Matrix<double> observation_operator(100, 100);
    observation_operator.Fill(0.);
    for(int i = 0; i < 100; i++)
        observation_operator(i, i) = 1.;
    observation_operator.Write("configuration/matrix.dat");

    END;

    return 0;

}
