// Copyright (C) 2008-2010 INRIA
// Author(s): Marc Fragu
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_CXX

#include "ReducedOrderExtendedKalmanFilter.hxx"

#include "seldon/computation/solver/SparseCholeskyFactorisation.cxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::ReducedOrderExtendedKalmanFilter(string configuration_file):
        configuration_file_(configuration_file)
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_,
                                     Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ReducedOrderExtendedKalmanFilter
                                     ::StaticMessage);


        /***************************
         * Reads the configuration *
         ***************************/


        Ops configuration(configuration_file_);
        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

        /*** Model ***/

        configuration.Set("model.configuration_file", "", configuration_file,
                          model_configuration_file_);

        /*** Observation manager ***/

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file,
                          observation_configuration_file_);

        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration",
                          option_display_["show_iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("display.show_time", option_display_["show_time"]);

        /*** Assimilation options ***/

        configuration.Set("data_assimilation.analyze_first_step",
                          analyze_first_step_);


        /*** Ouput saver ***/

        configuration.
            SetPrefix("reduced_order_extended_kalman_filter"
                      ".output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("state_forecast");
        output_saver_.Empty("state_analysis");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration",
                              output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::~ReducedOrderExtendedKalmanFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(bool initialize_model, bool initialize_observation_manager)
    {

        MessageHandler::Send(*this, "all", "::Initialize begin");

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);
        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

#ifdef VERDANDI_DENSE
        model_.GetStateErrorVarianceSqrt(L_, U_);
#else
        model_state_error_variance L_sparse, U_sparse;
        model_.GetStateErrorVarianceSqrt(L_sparse, U_sparse);
        Matrix<T, General, ArrayRowSparse> L_array, U_array;
        ConvertRowSparseToArrayRowSparse(L_sparse, L_array);
        ConvertRowSparseToArrayRowSparse(U_sparse, U_array);
        ConvertArrayRowSparseToDense(L_array, L_);
        ConvertArrayRowSparseToDense(U_array, U_);
#endif
        Nreduced_ = U_.GetN();

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        if (initialize_model)
        {
            MessageHandler::Send(*this, "model", "initial condition");
            MessageHandler::Send(*this, "driver", "initial condition");
        }

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the extended Kalman filter.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::InitializeStep()
    {
        model_.InitializeStep();
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        time_ = model_.GetTime();

        model_.Forward();

        PropagateCovarianceMatrix();

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Analyze()
    {
        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            if (option_display_["show_time"])
                cout << "Performing Reduced Order EKF at time step ["
                     << model_.GetTime() << "]..." << endl;

            model_state x;
            model_.GetState(x);
            Nstate_ = model_.GetNstate();

            observation y;
            observation_manager_.GetInnovation(x, y);
            Nobservation_ = y.GetSize();

            /*** Updated matrix U ***/

            dense_matrix HL(Nobservation_, Nreduced_),
                working_matrix_or(Nobservation_, Nreduced_);
            MltAdd(T(1), observation_manager_.GetTangentLinearOperator(),
                   L_, T(0), HL);
            MltAdd(T(1), observation_manager_.GetErrorVarianceInverse(), HL,
                   T(0), working_matrix_or);
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, working_matrix_or,
                   T(1), U_);

            /*** Updated K ***/

            dense_matrix U_inv(U_), K(Nstate_, Nobservation_),
                working_matrix_ro(Nreduced_, Nobservation_),
                working_matrix_ro2(Nreduced_, Nobservation_);

            GetInverse(U_inv);

#ifdef VERDANDI_DENSE
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans,
                   observation_manager_.GetErrorVarianceInverse(),
                   T(0), working_matrix_ro);
#else
            dense_matrix R_inv_dense;
            ConvertRowSparseToDense(
                observation_manager_.GetErrorVarianceInverse(), R_inv_dense);
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv_dense,
                   T(0), working_matrix_ro);
#endif
            MltAdd(T(1), U_inv, working_matrix_ro, T(0), working_matrix_ro2);
            MltAdd(T(1), L_, working_matrix_ro2, T(0), K);

            /*** Updated x ***/

            MltAdd(T(1), K, y, T(1), x);

            model_.SetState(x);

            if (option_display_["show_time"])
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Computes Covariance.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::PropagateCovarianceMatrix()
    {
        double saved_time = model_.GetTime();
        model_.SetTime(time_);

        // One column of L.
        model_state L_col(Nstate_);
        for (int j = 0; j < Nreduced_; j++)
        {
            GetCol(L_, j, L_col);
            model_.ApplyTangentLinearOperator(L_col);
            SetCol(L_col, j, L_);
        }

        model_.SetTime(saved_time);
    }



    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::HasFinished()
        const
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    const Model&
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::GetModel() const
    {
        return model_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::GetName() const
    {
        return "ReducedOrderExtendedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Message(string message)
    {
        model_state state;
        if (message.find("initial condition") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, double(model_.GetTime()),
                               "state_forecast");
        }

        if (message.find("forecast") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, double(model_.GetTime()),
                               "state_forecast");
        }

        if (message.find("analysis") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, double(model_.GetTime()),
                               "state_analysis");
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_CXX
#endif
