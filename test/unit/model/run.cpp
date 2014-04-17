#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT
#define VERDANDI_WITH_TRAJECTORY_MANAGER

#include <iostream>
#include <algorithm>
#include <cppunit/TestResult.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <tr1/random>

#include "Verdandi.hxx"
#include "seldon/computation/optimization/NLoptSolver.cxx"
#include "seldon/SeldonSolver.hxx"

#include "model/QuadraticModel.cxx"

using namespace CppUnit;
#ifdef VERDANDI_TEST_TIME
#include "test_time.hpp"
CPPUNIT_TEST_SUITE_REGISTRATION(TestTimeUnit);
#endif


int main(int argc, char** argv)
{
  VERDANDI_TRY;

  cout << endl ;
  TextUi::TestRunner runner;

  TestFactoryRegistry &registry = TestFactoryRegistry::getRegistry();

  runner.addTest(registry.makeTest());

  return runner.run("", false);

  VERDANDI_END;
}
