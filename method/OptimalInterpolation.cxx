// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_OPTIMALINTERPOLATION_CXX


#include "OptimalInterpolation.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::OptimalInterpolation(string configuration_file):
        model_(configuration_file),
        observation_manager_(model_, configuration_file),
        data_to_save_(false), analyzed_data_to_save_(false)
    {
        ConfigStream configuration_stream(configuration_file);

        /*** Initializations ***/

        model_.Initialize(configuration_file);
        data_to_save_ = true;
        observation_manager_.Initialize(model_, configuration_file);


        /***********************
         * Reads configuration *
         ***********************/


        /*** Display options ***/

        configuration_stream.SetSection("[display]");
        // Should iterations be displayed on screen?
        configuration_stream.PeekValue("Show_iterations",
                                       option_display_["show_iterations"]);
        // Should current date be displayed on screen?
        configuration_stream.PeekValue("Show_date",
                                       option_display_["show_date"]);

        /*** Assimilation options ***/

        Nstate_ = model_.GetNstate();

        configuration_stream.SetSection("[data_assimilation]");
        configuration_stream.PeekValue("Analyze_first_step",
                                       analyze_first_step_);

    }


    //! Destructor.
    template <class T, class ClassModel, class ClassObservationManager>
    OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::~OptimalInterpolation()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the optimal interpolation driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step.
      \param[in] configuration_file configuration file.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::Initialize(string configuration_file)
    {
        cout.precision(20);

        Vector<T> state_vector;

        /*** Initializations ***/

        // model_.Initialize(configuration_file);
        // observation_manager_.Initialize(model_, configuration_file);

        if (analyze_first_step_)
        {
            // Retrieves observations.
            observation_manager_.LoadObservation(model_);

            if (observation_manager_.IsAvailable())
            {
                if (option_display_["show_date"])
                    cout << "Performing optimal interpolation at time step ["
                         << model_.GetCurrentDate() << "]..." << endl;

                model_.GetState(state_vector);

                ComputeBLUE(state_vector);

                model_.SetState(state_vector);

                data_to_save_ = true;
                analyzed_data_to_save_ = true;

                if (option_display_["show_date"])
                    cout << " done." << endl;
            }
        }
    }


    //! Initializes a step for the optimal interpolation.
    /*! Initializes a step for the model.
     */
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::InitializeStep()
    {
        if (option_display_["show_date"])
            cout << "Current step: "
                 << model_.GetCurrentDate() << endl;
        model_.InitializeStep();
    }


    //! Performs a step forward without optimal interpolation.
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::Forward()
    {
        model_.Forward();
        data_to_save_ = true;
    }


    //! Computes the analyze.
    /*! To be called after the Forward method. Whenever observations are
      available, it assimilates them using optimal interpolation.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::Analyze()
    {
        Vector<T> state_vector;
        observation_manager_.LoadObservation(model_);

        if (observation_manager_.IsAvailable())
        {
            if (option_display_["show_date"])
                cout << "Performing optimal interpolation at time step ["
                     << model_.GetCurrentDate() << "]..." << endl;

            model_.GetState(state_vector);

            ComputeBLUE(state_vector);
            model_.SetState(state_vector);

            data_to_save_ = true;
            analyzed_data_to_save_ = true;

            if (option_display_["show_date"])
                cout << " done." << endl;
        }

    }


    //! Computes BLUE for optimal interpolation.
    /*! The state is updated by the combination of background state and
      innovation. It computes the BLUE (best linear unbiased estimator).
      \param[in] state_vector the state_vector to analyze.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::ComputeBLUE(Vector<T>& state_vector)
    {
        int r, c;

        // Number of observations at current date.
        Nobservation_ = observation_manager_.GetNobservation();

        // One row of background matrix B.
        Vector<T> error_covariance_row(Nstate_);

        // One row of tangent operator matrix.
        Vector<T> tangent_operator_row(Nstate_);

        // Temporary matrix and vector.
        Matrix<T> working_matrix(Nobservation_, Nobservation_);
        working_matrix.Fill(T(0));

        Vector<T> row(Nobservation_);

        // Computes HBH'.
        T H_entry;
        for (int j = 0; j < Nstate_; j++)
        {
            model_.GetBackgroundErrorCovarianceRow(j, error_covariance_row);
            // Computes the j-th row of BH'.
            for (r = 0; r < Nobservation_; r++)
            {
                observation_manager_
                    .GetTangentOperatorRow(r, tangent_operator_row);
                row(r) = DotProd(error_covariance_row, tangent_operator_row);
            }

            // Keeps on building HBH'.
            for (r = 0; r < Nobservation_; r++)
            {
                H_entry = observation_manager_.GetTangentOperator(r, j);
                for (c = 0; c < Nobservation_; c++)
                    working_matrix(r, c) += H_entry * row(c);
            }
        }

        // Computes (HBH' + R).
        for (r = 0; r < Nobservation_; r++)
            for (c = 0; c < Nobservation_; c++)
                working_matrix(r, c) += observation_manager_
                    .GetObservationErrorCovariance(r, c);

        // Computes (HBH' + R)^{-1}.
        GetInverse(working_matrix);

        // Computes innovation.
        Vector<T> innovation(Nobservation_);
        observation_manager_.GetInnovation(state_vector, innovation);

        // Computes working_matrix * innovation.
        Vector<T> working_vector(Nobservation_);
        MltAdd(T(1), working_matrix, innovation, T(0), working_vector);

        // Computes new state.
        Vector<T> working_vector0(Nobservation_);
        working_vector0.Fill(T(0));
        for (r = 0; r < Nstate_; r++)
        {
            // Computes the r-th row of BH'.
            model_.GetBackgroundErrorCovarianceRow(r, error_covariance_row);
            for (c = 0; c < Nobservation_; c++)
            {
                observation_manager_
                    .GetTangentOperatorRow(c, tangent_operator_row);
                working_vector0(c) = DotProd(error_covariance_row,
                                             tangent_operator_row);
            }

            state_vector(r) += DotProd(working_vector0, working_vector);
        }
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    bool OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::HasFinished() const
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    const ClassModel&
    OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::GetModel() const
    {
        return model_;
    }


    //! Checks if there is data to be saved.
    /*!
      \return True if there is data to be saved, false otherwise.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    bool OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::GetDataToSave() const
    {
        return data_to_save_;
    }


    //! Checks if there is analyzed data to be saved.
    /*!
      \return True if there is analyzed data to be saved, false otherwise.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    bool OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::GetAnalyzedDataToSave() const
    {
        return analyzed_data_to_save_;
    }


    //! Stores that there is no more data to be saved.
    template <class T, class ClassModel, class ClassObservationManager>
    void OptimalInterpolation<T, ClassModel, ClassObservationManager>
    ::ClearDataToSave()
    {
        data_to_save_ = false;
        analyzed_data_to_save_ = false;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OPTIMALINTERPOLATION_CXX
#endif
