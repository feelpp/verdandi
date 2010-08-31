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


#ifndef VERDANDI_FILE_METHOD_UNSCENTEDKALMANFILTER_CXX

#include "UnscentedKalmanFilter.hxx"

#include "BLUE.cxx"

#include "SigmaPoint.cxx"

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
    UnscentedKalmanFilter<T, Model, ObservationManager>
    ::UnscentedKalmanFilter(string configuration_file):
        model_(configuration_file),
        observation_manager_(model_, configuration_file)
    {
        Ops::Ops configuration(configuration_file);

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     UnscentedKalmanFilter::StaticMessage);


        /***************************
         * Reads the configuration *
         ***************************/


        configuration.SetPrefix("unscented_kalman_filter.");

        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration",
                          option_display_["show_iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("display.show_time", option_display_["show_time"]);

        /*** Assimilation options ***/

        configuration.Set("data_assimilation.analyze_first_step",
                          analyze_first_step_);

        /*** Sigma-points ***/

        configuration.Set("sigma_point.type",
                          "ops_in(v, {'canonical', 'star', 'simplex'})",
                          sigma_point_type_);

        /*** Ouput saver ***/

        configuration.SetPrefix("unscented_kalman_filter.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("state_forecast");
        output_saver_.Empty("state_analysis");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("unscented_kalman_filter.");

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
    template <class T, class Model, class ObservationManager>
    UnscentedKalmanFilter<T, Model, ObservationManager>
    ::~UnscentedKalmanFilter()
    {
        sigma_point_collection_.Deallocate();
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    void UnscentedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(string configuration_file)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        /*** Initializations ***/

        model_.Initialize(configuration_file);
        observation_manager_.Initialize(model_, configuration_file);

        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

        background_error_variance_.
            Copy(model_.GetStateErrorVariance());

        /*** Sigma-points ***/

        if (sigma_point_type_ == "canonical")
            ComputeCanonicalSigmaPoint(Nstate_, sigma_point_collection_,
                                       alpha_i_, alpha_constant_);
        else if (sigma_point_type_ == "star")
            ComputeStarSigmaPoint(Nstate_, sigma_point_collection_, alpha_i_,
                                  alpha_constant_);
        else if (sigma_point_type_ == "simplex")
            ComputeSimplexSigmaPoint(Nstate_, sigma_point_collection_,
                                     alpha_i_, alpha_constant_);

        if (alpha_constant_)
            alpha_ = alpha_i_(0);

        Nsigma_point_ = sigma_point_collection_.GetNvector();

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the unscented Kalman filter.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void UnscentedKalmanFilter<T, Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void UnscentedKalmanFilter<T, Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        // Computes of background error variance Cholesky factorization.
        model_state_error_variance background_error_variance_sqrt;
#ifdef VERDANDI_DENSE
        Copy(background_error_variance_, background_error_variance_sqrt);
        GetCholesky(background_error_variance_sqrt);
#else
        Matrix<T> tmp_sqrt;
        convert_RowSparse_to_Dense(background_error_variance_, tmp_sqrt);
        GetCholesky(tmp_sqrt);
        Matrix<T, General, ArrayRowSparse> tmp_sqrt_array;
        convert_Dense_to_ArrayRowSparse(tmp_sqrt, tmp_sqrt_array);
        Copy(tmp_sqrt_array, background_error_variance_sqrt);
#endif

        // Computes X_n^{(i)+}.
        model_state x;
        state_collection x_i;
        for (int i = 0; i < Nsigma_point_; i++)
        {
            model_.GetState(x);
            x_i.AddVector(x);
            x.Nullify();
            MltAdd(T(1), background_error_variance_sqrt,
                   sigma_point_collection_.GetVector(i), T(1),
                   x_i.GetVector(i));
        }


        if (alpha_constant_)
        {
            // Computes X_{n + 1}^-.
            x.Reallocate(Nstate_);
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                model_.ApplyOperator(x_i.GetVector(i), i + 1 == Nsigma_point_,
                                     true);
                Add(T(1), x_i.GetVector(i), x);
            }
            Mlt(alpha_, x);
            model_.SetState(x);

            // Computes P_{n + 1}^-.
#ifdef VERDANDI_DENSE
            background_error_variance_.Reallocate(Nstate_, Nstate_);
            background_error_variance_.Fill(T(0));
#else
            Matrix<T, General, ArrayRowSparse> tmp_array(Nstate_, Nstate_);
            Copy(tmp_array, background_error_variance_);
#endif
            model_state_error_variance working_matrix(Nstate_, 1);
            for (int i = 0; i < Nsigma_point_; i++)
            {
                Add(T(-1), x, x_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(x_i.GetVector(i), 0, working_matrix);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(x_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0, working_matrix);
#endif
                MltAdd(T(1), SeldonNoTrans, working_matrix, SeldonTrans,
                       working_matrix, T(1), background_error_variance_);
            }
            Mlt(alpha_, background_error_variance_);
        }
        else
        {
            // Computes X_{n + 1}^-.
            x.Reallocate(Nstate_);
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                model_.ApplyOperator(x_i.GetVector(i), i + 1 == Nsigma_point_,
                                     true);
                Add(alpha_i_(i), x_i.GetVector(i), x);
            }
            model_.SetState(x);

            // Computes P_{n + 1}^-.
#ifdef VERDANDI_DENSE
            background_error_variance_.Reallocate(Nstate_, Nstate_);
            background_error_variance_.Fill(T(0));
#else
            Matrix<T, General, ArrayRowSparse> tmp_array(Nstate_, Nstate_);
            Copy(tmp_array, background_error_variance_);
#endif
            model_state_error_variance working_matrix(Nstate_, 1);
            for (int i = 0; i < Nsigma_point_; i++)
            {
                Add(T(-1), x, x_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(x_i.GetVector(i), 0, working_matrix);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(x_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0, working_matrix);
#endif
                MltAdd(alpha_i_(i), SeldonNoTrans, working_matrix,
                       SeldonTrans, working_matrix, T(1),
                       background_error_variance_);
            }
        }

        x_i.Deallocate();

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE. */
    template <class T, class Model, class ObservationManager>
    void UnscentedKalmanFilter<T, Model, ObservationManager>::Analyze()
    {

        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (!observation_manager_.HasObservation())
        {
            MessageHandler::Send(*this, "all", "::Analyze end");
            return;
        }

        if (option_display_["show_time"])
            cout << "Performing UKF at time step ["
                 << model_.GetTime() << "]..." << endl;

        // Computes background error variance Cholesky factorization.
        model_state_error_variance background_error_variance_sqrt;
#ifdef VERDANDI_DENSE
        Copy(background_error_variance_, background_error_variance_sqrt);
        GetCholesky(background_error_variance_sqrt);
#else
        Matrix<T> tmp_sqrt;
        convert_RowSparse_to_Dense(background_error_variance_, tmp_sqrt);
        GetCholesky(tmp_sqrt);
        Matrix<T, General, ArrayRowSparse> tmp_sqrt_array;
        convert_Dense_to_ArrayRowSparse(tmp_sqrt, tmp_sqrt_array);
        Copy(tmp_sqrt_array, background_error_variance_sqrt);
#endif

        // Computes X_{n + 1}^{(i)-}.
        model_state x;
        state_collection x_i;
        for (int i = 0; i < Nsigma_point_; i++)
        {
            model_.GetState(x);
            x_i.AddVector(x);
            x.Nullify();
            MltAdd(T(1), background_error_variance_sqrt,
                   sigma_point_collection_.GetVector(i), T(1),
                   x_i.GetVector(i));
        }

        // Computes Z_{n + 1}^(i).
        Nobservation_ = observation_manager_.GetNobservation();
        model_state working_vector;
        state_collection z_i;
        for (int i = 0; i < Nsigma_point_; i++)
        {
            working_vector.Reallocate(Nobservation_);
            z_i.AddVector(working_vector);
            working_vector.Nullify();
            observation_manager_.ApplyOperator(x_i.GetVector(i),
                                               z_i.GetVector(i));
        }

        if (alpha_constant_)
        {
            // Computes the predicted measurement Z_{n + 1}.
            model_state z;
            z.Reallocate(Nobservation_);
            z.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
                Add(T(1), z_i.GetVector(i), z);
            Mlt(alpha_, z);

            // Computes X_{n+1}-.
            x.Reallocate(Nstate_);
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
                Add(T(1), x_i.GetVector(i), x);
            Mlt(alpha_, x);

            // Computes P_XZ = cov(X_{n + 1}^*, Z_{n + 1}^*).
            model_state_error_variance P_xz(Nstate_, Nobservation_);
#ifdef VERDANDI_DENSE
            P_xz.Fill(T(0));
#else
            Matrix<T, General, ArrayRowSparse> tmp_array(Nstate_, Nstate_);
            Copy(tmp_array, P_xz);
#endif
            model_state_error_variance working_matrix_state(Nstate_, 1),
                working_matrix_observation(Nobservation_, 1);
            for (int i = 0; i < Nsigma_point_; i++)
            {
                Add(T(-1), x, x_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(x_i.GetVector(i), 0, working_matrix_state);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(x_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0, working_matrix_state);
#endif
                Add(T(-1), z, z_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(z_i.GetVector(i), 0, working_matrix_observation);
#else
                model_state_error_variance_row
                    error_covariance_column2(Nstate_);
                convert_Dense_to_Sparse(z_i.GetVector(i),
                                        error_covariance_column2);
                SetCol(error_covariance_column2, 0,
                       working_matrix_observation);
#endif
                MltAdd(T(1), SeldonNoTrans, working_matrix_state,
                       SeldonTrans, working_matrix_observation,
                       T(1), P_xz);
            }
            Mlt(alpha_, P_xz);

            // Computes P_Z = cov(Z_{n + 1}^*, Z_{n + 1}^*).
            model_state_error_variance P_z(Nobservation_, Nobservation_);
            P_z.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
#ifdef VERDANDI_DENSE
                SetCol(z_i.GetVector(i), 0, working_matrix_observation);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(z_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0,
                       working_matrix_observation);
#endif
                MltAdd(T(1), SeldonNoTrans,
                       working_matrix_observation, SeldonTrans,
                       working_matrix_observation, T(1), P_z);
            }
            Mlt(alpha_, P_z);
            Add(T(1), observation_manager_.GetErrorVariance(), P_z);

            // Computes the Kalman gain K_{n + 1}.
            matrix_state_observation K(Nstate_, Nobservation_);
            K.Fill(T(0));
            GetInverse(P_z);
            MltAdd(T(1), P_xz, P_z, T(0), K);

            // Computes X_{n + 1}^+.
            model_state x0;
            model_.GetState(x0);
            observation innovation;
            observation_manager_.GetInnovation(x, innovation);
            MltAdd(T(1), K, innovation, T(1), x0);
            model_.SetState(x0);

            // Computes P_{n + 1}^+.
            MltAdd(T(-1), SeldonNoTrans, K, SeldonTrans, P_xz, T(1),
                   background_error_variance_);
        }
        else
        {
            // Computes the predicted measurement Z_{n + 1}.
            model_state z;
            z.Reallocate(Nobservation_);
            z.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
                Add(alpha_i_(i), z_i.GetVector(i), z);

            // Computes X_{n+1}-.
            x.Reallocate(Nstate_);
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
                Add(alpha_i_(i), x_i.GetVector(i), x);

            // Computes P_XZ = cov(X_{n + 1}^*, Z_{n + 1}^*).
            model_state_error_variance P_xz(Nstate_, Nobservation_);
#ifdef VERDANDI_DENSE
            P_xz.Fill(T(0));
#else
            Matrix<T, General, ArrayRowSparse> tmp_array(Nstate_, Nstate_);
            Copy(tmp_array, P_xz);
#endif
            model_state_error_variance working_matrix_state(Nstate_, 1),
                working_matrix_observation(Nobservation_, 1);
            for (int i = 0; i < Nsigma_point_; i++)
            {
                Add(T(-1), x, x_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(x_i.GetVector(i), 0, working_matrix_state);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(x_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0, working_matrix_state);
#endif
                Add(T(-1), z, z_i.GetVector(i));
#ifdef VERDANDI_DENSE
                SetCol(z_i.GetVector(i), 0, working_matrix_observation);
#else
                model_state_error_variance_row
                    error_covariance_column2(Nstate_);
                convert_Dense_to_Sparse(z_i.GetVector(i),
                                        error_covariance_column2);
                SetCol(error_covariance_column2, 0,
                       working_matrix_observation);
#endif
                MltAdd(alpha_i_(i), SeldonNoTrans, working_matrix_state,
                       SeldonTrans, working_matrix_observation, T(1), P_xz);
            }

            // Computes P_Z = cov(Z_{n + 1}^*, Z_{n + 1}^*).
            model_state_error_variance P_z(Nobservation_, Nobservation_);
            P_z.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
#ifdef VERDANDI_DENSE
                SetCol(z_i.GetVector(i), 0, working_matrix_observation);
#else
                model_state_error_variance_row
                    error_covariance_column(Nstate_);
                convert_Dense_to_Sparse(z_i.GetVector(i),
                                        error_covariance_column);
                SetCol(error_covariance_column, 0,
                       working_matrix_observation);
#endif
                MltAdd(alpha_i_(i), SeldonNoTrans,
                       working_matrix_observation, SeldonTrans,
                       working_matrix_observation, T(1), P_z);
            }
            Add(T(1), observation_manager_.GetErrorVariance(), P_z);

            // Computes the Kalman gain K_{n + 1}.
            matrix_state_observation K(Nstate_, Nobservation_);
            K.Fill(T(0));
            GetInverse(P_z);
            MltAdd(T(1), P_xz, P_z, T(0), K);

            // Computes X_{n + 1}^+.
            model_state x0;
            model_.GetState(x0);
            observation innovation;
            observation_manager_.GetInnovation(x, innovation);
            MltAdd(T(1), K, innovation, T(1), x0);
            model_.SetState(x0);

            // Computes P_{n + 1}^+.
            MltAdd(T(-1), SeldonNoTrans, K, SeldonTrans, P_xz, T(1),
                   background_error_variance_);
        }

        x_i.Deallocate();
        z_i.Deallocate();

        if (option_display_["show_time"])
            cout << " done." << endl;

        MessageHandler::Send(*this, "model", "analysis");
        MessageHandler::Send(*this, "observation_manager", "analysis");
        MessageHandler::Send(*this, "driver", "analysis");

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool UnscentedKalmanFilter<T, Model, ObservationManager>::HasFinished()
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
    UnscentedKalmanFilter<T, Model, ObservationManager>::GetModel() const
    {
        return model_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string
    UnscentedKalmanFilter<T, Model, ObservationManager>::GetName() const
    {
        return "UnscentedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void UnscentedKalmanFilter<T, Model, ObservationManager>
    ::Message(string message)
    {
        model_state state;
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


#define VERDANDI_FILE_METHOD_UNSCENTEDKALMANFILTER_CXX
#endif
