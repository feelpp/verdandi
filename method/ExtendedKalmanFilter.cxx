// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef VERDANDI_FILE_EXTENDEDKALMANFILTER_CXX


#include "ExtendedKalmanFilter.hxx"

#include "BLUE.cxx"

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
    ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::ExtendedKalmanFilter(string configuration_file):
        model_(configuration_file),
        observation_manager_(model_, configuration_file)
    {
        Ops::Ops configuration(configuration_file);

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_,
                                     ClassModel::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ClassObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ExtendedKalmanFilter::StaticMessage);


        /***************************
         * Reads the configuration *
         ***************************/


        configuration.SetPrefix("extended_kalman_filter.");

        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration",
                          option_display_["show_iteration"]);
        // Should current date be displayed on screen?
        configuration.Set("display.show_date", option_display_["show_date"]);

        /*** Assimilation options ***/

        configuration.Set("data_assimilation.analyze_first_step",
                          analyze_first_step_);

        configuration.Set("BLUE_computation",
                          "ops_in(v, {'vector', 'matrix'})",
                          blue_computation_);
        configuration.Set("covariance_computation",
                          "ops_in(v, {'vector', 'matrix'})",
                          covariance_computation_);

        /*** Ouput saver ***/

        configuration.SetPrefix("extended_kalman_filter.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("state_forecast");
        output_saver_.Empty("state_analysis");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("extended_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }
    }


    //! Destructor.
    template <class T, class ClassModel, class ClassObservationManager>
    ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::~ExtendedKalmanFilter()
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
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::Initialize(string configuration_file)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        /*** Initializations ***/

        model_.Initialize(configuration_file);
        observation_manager_.Initialize(model_, configuration_file);

        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

        background_error_variance_.
            Copy(model_.GetBackgroundErrorVarianceMatrix());

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the optimal interpolation.
    /*! Initializes a step for the model.
     */
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        model_.Forward();

        PropagateCovarianceMatrix();

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE.
     */
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::Analyze()
    {

        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetDate(model_, model_.GetDate());

        if (observation_manager_.HasObservation())
        {
            if (option_display_["show_date"])
                cout << "Performing EKF at time step ["
                     << model_.GetDate() << "]..." << endl;

            state_vector state;
            model_.GetState(state);
            Nstate_ = model_.GetNstate();

            observation_vector innovation;
            observation_manager_.GetInnovation(state, innovation);
            Nobservation_ = innovation.GetSize();

            ComputeBLUE(innovation, state);

            model_.SetState(state);

            if (option_display_["show_date"])
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Computes covariance.
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::PropagateCovarianceMatrix_vector()
    {
        double saved_date;
        state_vector saved_state;
        saved_date = model_.GetDate();
        model_.GetState(saved_state);

        // One column of covariance matrix P.
        background_error_covariance_vector error_covariance_column(Nstate_);

        for (int j = 0; j < Nstate_; j++)
        {
            GetCol(background_error_variance_, j,
                   error_covariance_column);
#ifdef VERDANDI_DENSE
            model_.ApplyTangentLinearModel(error_covariance_column);
#else
            state_vector working_vector(Nstate_);
            ConvertSparseToDense(error_covariance_column, working_vector);
            model_.ApplyTangentLinearModel(working_vector);
            ConvertDenseToSparse(working_vector, error_covariance_column);
#endif
            SetCol(error_covariance_column, j, background_error_variance_);
        }

        Transpose(background_error_variance_);

        for (int j = 0; j < Nstate_; j++)
        {
            GetCol(background_error_variance_, j,
                   error_covariance_column);
#ifdef VERDANDI_DENSE
            model_.ApplyTangentLinearModel(error_covariance_column);
#else
            state_vector working_vector(Nstate_);
            ConvertSparseToDense(error_covariance_column, working_vector);
            model_.ApplyTangentLinearModel(working_vector);
            ConvertDenseToSparse(working_vector, error_covariance_column);
#endif
            SetCol(error_covariance_column, j, background_error_variance_);
        }

        model_.SetDate(saved_date);
        model_.SetState(saved_state);
    }


    //! Computes covariance.
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::PropagateCovarianceMatrix_matrix()
    {
        tangent_operator_matrix A;
        model_.GetTangentLinearModel(A);

        MltAdd(T(1.), A, background_error_variance_,
               T(0.), background_error_variance_);

        MltAdd(T(1.), SeldonNoTrans,
               background_error_variance_, SeldonTrans, A, T(0.),
               background_error_variance_);
    }


    //! Computes BLUE for Extended Kalman Filter.
    /*! The state is updated by the combination of background state and
      innovation. It computes the BLUE (best linear unbiased estimator).
      \param[in] state_vector the state vector to analyze.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::ComputeBLUE(const observation_vector& innovation, state_vector& state)
    {
        if (blue_computation_ == "vector")
            throw ErrorUndefined("ExtendedKalmanFilter"
                                 "::ComputeBLUE()");
        else
            ComputeBLUE_matrix(background_error_variance_,
                               observation_manager_.
                               GetTangentOperatorMatrix(),
                               innovation,
                               observation_manager_.
                               GetObservationErrorVariance(),
                               state, true, true);
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    bool ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::HasFinished() const
    {
        return model_.HasFinished();
    }


    //! Computes Covariance.
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::PropagateCovarianceMatrix()
    {
        if (covariance_computation_ == "vector")
            PropagateCovarianceMatrix_vector();
        else
            PropagateCovarianceMatrix_matrix();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    const ClassModel&
    ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::GetModel() const
    {
        return model_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    string ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::GetName() const
    {
        return "ExtendedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class ClassModel, class ClassObservationManager>
    void ExtendedKalmanFilter<T, ClassModel, ClassObservationManager>
    ::Message(string message)
    {
        state_vector state;
        if (message.find("forecast") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, double(model_.GetDate()),
                               "state_forecast");
        }

        if (message.find("analysis") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, double(model_.GetDate()),
                               "state_analysis");
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_EXTENDEDKALMANFILTER_CXX
#endif
