// Copyright (C) 2014 INRIA
// Author(s): Nicolas Claude
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#include <cppunit/extensions/HelperMacros.h>
#include <iostream>

using namespace std;


class TestTimeUnit: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestTimeUnit);
    CPPUNIT_TEST(TestTimeBackward);
    CPPUNIT_TEST(TestTimeRandom);
    CPPUNIT_TEST_SUITE_END();
    typedef double real;
    typedef Verdandi::Vector<double> state;

protected:
    // The accuracy of the test.
    double accuracy_;
    // The size of the state.
    unsigned int Nstate_;
    // Generates random numbers, using Mersenne Twister.
    tr1::mt19937 generator_;
    // Tested model.
    VERDANDI_CPPUNIT_MODEL model_;
    // Stores the states to check the consistency of the model.
    vector<state> state_storage_;
    // Stores the times to check the consistency of the model.
    vector<double> time_storage_;
    // Span of the storage.
    int Nstep_;

public:
    // This method initializes the model, reads the configuration for the
    // tests and stores all the states and times needed for the tests.
    void setUp()
    {
        // Initializes the random generator.
        generator_.seed(time(NULL));
        // Intializes the model, reads the configuration, initializes
        // variables.
        model_.Initialize(VERDANDI_CPPUNIT_MODEL_CONFIG);
        Verdandi::VerdandiOps configuration(VERDANDI_CPPUNIT_MODEL_CONFIG);
        configuration.Set("Nstep_time", Nstep_);
        Nstate_ = model_.GetNstate();

        // Computes and stores the states and times.
        for (int i = 0; i < Nstep_; i++)
        {
            state_storage_.push_back(model_.GetState());
            time_storage_.push_back(model_.GetTime());
            model_.Forward();
        }
    }


    void tearDown()
    {
    }


    // This test checks that the model will give the correct state, given the
    // previous state.
    void TestTimeBackward()
    {
        cout << "Test time backward: ";
        // For each time step, asks the model to go forward in time once and
        // checks if the given state is correct.
        for (int i = 0; i < Nstep_ - 1; i++)
        {
            model_.GetState() = state_storage_[Nstep_ - i - 2];
            model_.StateUpdated();
            model_.SetTime(time_storage_[Nstep_ - i - 2]);
            model_.Forward();
            for(unsigned int j = 0; j < Nstate_; j++)
                CPPUNIT_ASSERT(model_.GetState()(j) ==
                               state_storage_[Nstep_ - i - 1](j));
        }
        cout << " OK" << endl << endl;
    }


    // This test takes two random points in time and tries to pull the model
    // from the first to the second.
    void TestTimeRandom()
    {
        cout << "Test time random: ";
        // Randomizes two steps.
        unsigned int step1 = 0;
        unsigned int step2 = 0;
        while (step1 == step2 || step2 < step1)
        {
            step1 = generator_() % Nstep_;
            step2 = generator_() % Nstep_;
        }

        // Inserts the state at starting time into the model.
        model_.GetState() = state_storage_[step1];
        model_.StateUpdated();

        // Inserts the correct time into the model.
        model_.SetTime(time_storage_[step1]);

        // Pulls the model forward in time.
        for(unsigned int i = step1; i < step2; i++)
            model_.Forward();

        // Checks if the computed state is correct.
        for(unsigned int j = 0; j < Nstate_; j++)
            CPPUNIT_ASSERT_MESSAGE("Discrepancy detected",
                                   model_.GetState()(j)
                                   == state_storage_[step2](j));

        cout << " OK" << endl << endl;
    }
};
