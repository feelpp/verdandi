// Copyright (C) 2008-2010 INRIA
// Author(s): Vivien Mallet, Serhiy Zhuk
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDMINIMAX_CXX


#include "ReducedMinimax.hxx"


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
    ReducedMinimax<T, Model, ObservationManager>
    ::ReducedMinimax(string configuration_file):
        configuration_file_(configuration_file)
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_,
                                     Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ReducedMinimax::StaticMessage);
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedMinimax<T, Model, ObservationManager>::~ReducedMinimax()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Initialize()
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        Ops configuration(configuration_file_);
        configuration.SetPrefix("reduced_minimax.");


        /*********************************
         * Model and observation manager *
         *********************************/


        configuration.Set("model.configuration_file",
                          model_configuration_file_);
        model_.Initialize(model_configuration_file_);
        Nstate_ = model_.GetNstate();

        configuration.Set("observation_manager.configuration_file",
                          observation_configuration_file_);
        observation_manager_.Initialize(model_,
                                        observation_configuration_file_);


        /***************************
         * Reads the configuration *
         ***************************/


        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration", show_iteration_);
        // Should current time be displayed on screen?
        configuration.Set("display.show_time", show_time_);

        iteration_ = 0;

        /*** Reduction options ***/

        configuration.Set("reduction_method", "ops_in(v, {'none', 'pod'})",
                          reduction_method_);
        if (reduction_method_ == "pod")
        {
            configuration.Set("pod.Nprojection_max", "v > 0",
                              Nprojection_max_);
            configuration.Set("pod.acceptable_error",
                              "v >= 0 and v <= 1",
                              acceptable_error_);
            configuration.Set("pod.Nsnapshot", "v > 1", Nsnapshot_);
            // There cannot be more projection modes than snapshots.
            Nprojection_max_ = min(Nprojection_max_, Nsnapshot_);
        }

        /*** Model options ***/

        configuration.Set("model.bound_over_standard_deviation",
                          bound_over_standard_deviation_);

        /*** Filter options ***/

        T diagonal_part;
        if (configuration.Is<T>("model_error.diagonal_part"))
        {
            configuration.Set("model_error.diagonal_part", diagonal_part);
            D_tilde_inv_.Reallocate(Nstate_);
            D_tilde_inv_.Fill(diagonal_part);
        }
        else
            configuration.Set("model_error.diagonal_part", D_tilde_inv_);
        // 'D_tilde_inv_' is provided as a variance.
        for (int i = 0; i < D_tilde_inv_.GetLength(); i++)
            D_tilde_inv_(i) = bound_over_standard_deviation_
                * sqrt(D_tilde_inv_(i));
        // 'D_tilde_inv_' is inverted.
        for (int i = 0; i < D_tilde_inv_.GetLength(); i++)
            D_tilde_inv_(i) = 1. / D_tilde_inv_(i);

        // Systematic error in the initial condition.
        T error;
        if (configuration.Is<T>("systematic_error.initial_condition"))
        {
            configuration.Set("systematic_error.initial_condition", error);
            e0_.Reallocate(Nstate_);
            e0_.Fill(error);
        }
        else
            configuration.Set("systematic_error.initial_condition", e0_);

        // Systematic error in the model error.
        e_.Reallocate(Nstate_);
        e_.Fill(0.);

        // Systematic error in the observation error.
        configuration.Set("systematic_error.observation_error", eta_);

        /*** Ouput saver ***/

        configuration.SetPrefix("reduced_minimax.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("full_analysis");
        output_saver_.Empty("reduced_analysis");
        output_saver_.Empty("projected_state");
        output_saver_.Empty("projection");
        output_saver_.Empty("minimax_gain");
        output_saver_.Empty("snapshot");
        output_saver_.Empty("singular_value");
        output_saver_.Empty("left_singular_vector");
        output_saver_.Empty("right_singular_vector");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("reduced_minimax.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }


        /*************
         * Algorithm *
         *************/

        if (reduction_method_ == "pod")
            // Starting with in a POD sequence.
            mode_ = 0;
        else
        {
            Nprojection_ = Nstate_;
            projection_.Reallocate(Nprojection_, Nprojection_);
            projection_.SetIdentity();
            previous_projection_ = projection_;
            Nprevious_projection_ = previous_projection_.GetM();
            mode_ = 1;
        }
        inner_iteration_ = 0;
        first_sequence_ = true;

        /*** Filter ***/

        if (show_time_)
            Logger::StdOut(*this,
                           "Starting time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Starting time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Initialization");
        else
            Logger::Log<-3>(*this, "Initialization");

        MessageHandler::Send(*this, "model", "initial condition");

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step.
    /*! Initializes a step for the model. */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        if (mode_ == 1 && inner_iteration_ == 0)
        {
            if (reduction_method_ != "none")
            {
                model_.SetFullState(full_state_);
                model_.SetTime(starting_time_);
            }

            if (show_time_ || show_iteration_)
                Logger::StdOut(*this, "Starting filtering sequence");
            else
                Logger::Log<-3>(*this, "Starting filtering sequence");

            /*** Filter ***/

            if (first_sequence_)
            {
                FilterInitialization();

                MessageHandler::Send(*this, "model", "full_analysis");
                MessageHandler::Send(*this, "observation_manager",
                                     "full_analysis");
                MessageHandler::Send(*this, "driver", "full_analysis");
                MessageHandler::Send(*this, "driver", "reduced_analysis");
            }

            // In other windows, there is no such special case. The error
            // associated with the first step in the window is the error
            // of the last step in the previous window.
            first_sequence_ = false;
        }

        model_.InitializeStep();

        if (mode_ == 0 && inner_iteration_ == 0)
        {
            // Saves the initial state for the next sequence of error
            // computation.
            starting_time_ = model_.GetTime();
            model_.GetFullState(full_state_);

            model_state state;
            model_.GetState(state);

            snapshot_.Reallocate(Nstate_, Nsnapshot_);
            SetCol(state, inner_iteration_++, snapshot_);

            if (show_time_ || show_iteration_)
                Logger::StdOut(*this, "Starting POD sequence");
            else
                Logger::Log<-3>(*this, "Starting POD sequence");
        }

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        if (show_time_)
            Logger::StdOut(*this, "Time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Iteration " + to_str(iteration_) + " -> "
                           + to_str(iteration_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(iteration_) + " -> "
                            + to_str(iteration_ + 1));

        if (mode_ == 0)
            SnapshotRecording();
        else
        {
            Propagation();

            MessageHandler::Send(*this, "model", "full_analysis");
            MessageHandler::Send(*this, "observation_manager",
                                 "full_analysis");
            MessageHandler::Send(*this, "driver", "full_analysis");
            MessageHandler::Send(*this, "driver", "reduced_analysis");
        }

        iteration_++;

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Initialization of the filter.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::FilterInitialization()
    {
        MessageHandler::Send(*this, "all", "::FilterInitialization begin");

        // Model state.
        double time = model_.GetTime();
        model_state full_state;
        model_.GetFullState(full_state);

        // Observation operator and data.
	observation_tangent_linear_operator H_tilde;
	observation_tangent_linear_operator H;
	Vector<T> y;

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            observation_manager_.GetObservation(y);
            Nobservation_ = y.GetLength();
            for (int i = 0; i < Nobservation_; i++)
                y(i) -= eta_;
            H_tilde = observation_manager_.GetTangentLinearOperator();
            H.Reallocate(Nobservation_, Nprojection_);
            MltAdd(T(1), SeldonNoTrans, H_tilde, SeldonTrans, projection_,
                   T(0), H);
            R_inv_ = observation_manager_.GetErrorVariance();
            GetInverse(R_inv_);
            Mlt(1. / (bound_over_standard_deviation_
                      * bound_over_standard_deviation_), R_inv_);
        }
        else
        {
            // If there are no observations, we add a fictitious observation
            // and a null H.
            y.Reallocate(1);
            y.Zero();
            Nobservation_ = 1;
            H.Reallocate(1, Nprojection_);
            H.Zero();
            R_inv_.Reallocate(1, 1);
            R_inv_(0, 0) = T(1);
        }

        // Temporary variables.
	Matrix<T, General, RowMajor> mtmp, mtmp_1;
        Vector<T> vtmp;

        /*** Computes the minimax gain 'G_' ***/

        // Computes $\widecheck F_0$.
	Matrix<T, General, RowMajor> F_check(Nstate_, Nprojection_);
        vtmp.Reallocate(Nprojection_);
        for (int i = 0; i < Nstate_; i++)
        {
            GetCol(projection_, i, vtmp);
            Mlt(D_tilde_inv_(i), vtmp);
            SetRow(vtmp, i, F_check);
        }

        // Computes $\widecheck Q^{\frac 12}$.
	model_state_error_variance Q_sqrt_check;
        Q_sqrt_check = model_.GetStateErrorVarianceSqrt();
        Nmode_Q_= Q_sqrt_check.GetN();
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nmode_Q_; j++)
                Q_sqrt_check(i, j) *= bound_over_standard_deviation_
                    * D_tilde_inv_(i);

        // Puts $\widecheck F_0^T \widecheck F_0$ into 'G_'.
        G_.Reallocate(Nprojection_, Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, F_check, T(0), G_);

        // Computes $\widecheck F_0^T \widecheck Q_^{\frac 12}$.
        Matrix<T, General, RowMajor> FtQ(Nprojection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, F_check,
               SeldonNoTrans, Q_sqrt_check, T(0), FtQ);

        // Computes $(I_{q\times q} + \widecheck Q^{\frac T2} \widecheck
        // Q^{\frac 12})^{-1}$.
	Matrix<T, General, RowMajor> IQtQinv(Nmode_Q_, Nmode_Q_);
        IQtQinv.SetIdentity();
        MltAdd(T(1), SeldonTrans, Q_sqrt_check,
               SeldonNoTrans, Q_sqrt_check, T(1), IQtQinv);
        // Note that IQtQinv is used for estimate's propagation too.
        GetInverse(IQtQinv);

        // Subtracts to 'G_' this part: $\widecheck F_{0}^T\widecheck Q^{\frac
        // 12} (I_{q\times q} + \widecheck Q^{\frac T2} \widecheck Q^{\frac
        // 12})^{-1} (\widecheck F_{0}^T \widecheck Q^{\frac 12})^T$.
        Matrix<T, General, RowMajor> FtQ_IQtQinv(Nprojection_, Nmode_Q_);
        MltAdd(T(1), FtQ, IQtQinv, T(0), FtQ_IQtQinv);
        MltAdd(T(-1), SeldonNoTrans, FtQ_IQtQinv, SeldonTrans, FtQ, T(1), G_);

        // Adds $H_0^T R_0^{-1} H_0$ to 'G_'.
        Matrix<T, General, RowMajor> Ht_Rinv(Nprojection_, Nobservation_);
        MltAdd(T(1), SeldonTrans, H, SeldonNoTrans, R_inv_, T(0), Ht_Rinv);
        MltAdd(T(1), Ht_Rinv, H, T(1), G_);

        /*** Computes the minimax estimator ***/

        Vector<T> e0;
        model_.GetState(e0);
        Add(T(1), e0_, e0);
        for (int j = 0; j < Nstate_; j++)
            e0(j) *= D_tilde_inv_(j);

        // Computes $F_0^T \widetilde D^{-\frac 12}\overline e$.
        z_.Reallocate(Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, e0, T(0), z_);

        // Computes $F_0^T \widecheck Q^{\frac 12} (I_{q \times q} +
        // \widecheck Q^{\frac T2} \widecheck Q^{\frac 12})^{-1} \widecheck
        // Q^{\frac T2} \widetilde D^{-\frac 12 } \overline e$.
        vtmp.Reallocate(Nmode_Q_);
        MltAdd(T(1), SeldonTrans, Q_sqrt_check, e0, T(0), vtmp);
        MltAdd(T(-1), FtQ_IQtQinv, vtmp, T(1), z_);

        // Adds $H^T_0 R_0^{-1} (y_0 - \eta_0)$.
        MltAdd(T(1), Ht_Rinv, y, T(1), z_);

        // Finally computes the minimax state, with $G_0^{-1} z_0$.
        state_ = z_;
        mtmp = G_;
        GetAndSolveLU(mtmp, state_);

        vtmp.Reallocate(Nstate_);
        MltAdd(T(1), SeldonTrans, projection_, state_, T(0), vtmp);
        model_.SetState(vtmp);

        output_saver_.Save(G_, model_.GetTime(), "minimax_gain");
        output_saver_.Save(projection_, "projection");

        MessageHandler::Send(*this, "all", "::FilterInitialization end");
    }


    //! Propagates the state and the minimax gain.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Propagation()
    {
        MessageHandler::Send(*this, "all", "::Propagation begin");

        // Model state.
        double time = model_.GetTime();
        model_state full_state;
        model_.GetFullState(full_state);

        // Tangent linear model.
	Matrix<T, General, RowMajor> M(Nstate_, Nprevious_projection_);
	ComputeTangentLinearModel(M);

        // Temporary variables.
        Vector<T> vtmp, vtmp_1;
	Matrix<T, General, RowMajor> mtmp, mtmp_1;

        /*** Model-related variables ***/

        // Computes $\widecheck M_t$.
	Matrix<T, General, RowMajor> M_check = M;
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nprevious_projection_; j++)
                M_check(i, j) *= D_tilde_inv_(i);

        // Model error.
        model_error_variance Q_sqrt;
        Q_sqrt = model_.GetErrorVarianceSqrt();
        Nmode_Q_= Q_sqrt.GetN();

        // Computes $\widecheck Q_t^{\frac 12}$.
	Matrix<T, General, RowMajor> Q_sqrt_check = Q_sqrt;
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nmode_Q_; j++)
                Q_sqrt_check(i, j) *= bound_over_standard_deviation_
                    * D_tilde_inv_(i);
        Q_sqrt.Clear();

        // Computes $G_t + \widecheck M_t^T \widecheck M_t$.
	Matrix<T, General, RowMajor> G_MtM(Nprevious_projection_,
                                           Nprevious_projection_);
        MltAdd(T(1), SeldonTrans, M_check,
               SeldonNoTrans, M_check, T(0), G_MtM);
        Add(T(1), G_, G_MtM);

        // Computes $\widecheck U_t$.
	Matrix<T, General, RowMajor> U_check(Nstate_, Nprevious_projection_);
        mtmp = G_MtM;
        GetInverse(mtmp);
        GetCholesky(mtmp);
        MltAdd(T(1), M_check, mtmp, T(0), U_check);

        // Computes $V_t$.
	Matrix<T, General, RowMajor> IQtQinv;
	Matrix<T, General, RowMajor> Vinv(Nmode_Q_, Nmode_Q_);
        Vinv.SetIdentity();
        MltAdd(T(1), SeldonTrans, Q_sqrt_check,
               SeldonNoTrans, Q_sqrt_check, T(1), Vinv);
        IQtQinv = Vinv;
        // Note that IQtQinv is used for estimate's propagation too.
        GetInverse(IQtQinv);

        mtmp.Reallocate(Nmode_Q_, Nprevious_projection_);
        MltAdd(T(1), SeldonTrans,
               Q_sqrt_check, SeldonNoTrans, U_check, T(0), mtmp);
        MltAdd(T(-1), SeldonNoTrans, mtmp, SeldonTrans, mtmp, T(1), Vinv);
        GetInverse(Vinv);

        /*** Computes $B_t$ and lets the model perform one time step ***/

        // Computes $B_t$.
        B_ = G_MtM;
        mtmp.Reallocate(Nprevious_projection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, M_check,
               SeldonNoTrans, Q_sqrt_check, T(0), mtmp);

        mtmp_1.Reallocate(Nprevious_projection_, Nmode_Q_);
        MltAdd(T(1), mtmp, IQtQinv, T(0), mtmp_1);

        MltAdd(T(-1), SeldonNoTrans, mtmp_1, SeldonTrans, mtmp, T(1), B_);
        GetInverse(B_);

        // Calls the model to compute $\mathcal{\widetilde M}_t(F_{t+1} B_t
        // z_t)$.
        vtmp.Reallocate(Nprevious_projection_);
        Mlt(B_, z_, vtmp);
        Vector<T> MFBz(Nstate_);
        MltAdd(T(1), SeldonTrans, previous_projection_, vtmp, T(0), MFBz);
        model_.ApplyOperator(MFBz, true, true);
        // Note that the model is now one time step further.

        /*** Observation-related variables ***/

        // Observation operator and data.
	observation_tangent_linear_operator H_tilde;
	observation_tangent_linear_operator H;
	Vector<T> y;

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            observation_manager_.GetObservation(y);
            Nobservation_ = y.GetLength();
            for (int i = 0; i < Nobservation_; i++)
                y(i) -= eta_;
            H_tilde = observation_manager_.GetTangentLinearOperator();
            H.Reallocate(Nobservation_, Nprojection_);
            MltAdd(T(1), SeldonNoTrans, H_tilde, SeldonTrans, projection_,
                   T(0), H);
            R_inv_ = observation_manager_.GetErrorVariance();
            GetInverse(R_inv_);
            Mlt(1. / (bound_over_standard_deviation_
                      * bound_over_standard_deviation_), R_inv_);
        }
        else
        {
            // If there are no observations, we add a fictitious observation
            // and a null H.
            y.Reallocate(1);
            y.Zero();
            Nobservation_ = 1;
            H.Reallocate(1, Nprojection_);
            H.Zero();
            R_inv_.Reallocate(1, 1);
            R_inv_(0, 0) = T(1);
        }

        /*** Computes the minimax gain $G_t$ ***/

        // Computes $\widecheck F_{t+1}$.
	Matrix<T, General, RowMajor> F_check(Nstate_, Nprojection_);
        vtmp.Reallocate(Nprojection_);
        for (int i = 0; i < Nstate_; i++)
        {
            GetCol(projection_, i, vtmp);
            Mlt(D_tilde_inv_(i), vtmp);
            SetRow(vtmp, i, F_check);
        }

        // Computes $G_t$.
        G_.Reallocate(Nprojection_, Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, F_check, T(0), G_);

	Matrix<T, General, RowMajor> FtU(Nprojection_, Nprevious_projection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, U_check, T(0), FtU);
        MltAdd(T(-1), SeldonNoTrans, FtU, SeldonTrans, FtU, T(1), G_);

	Matrix<T, General, RowMajor> FtQ(Nprojection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, F_check,
               SeldonNoTrans, Q_sqrt_check, T(0), FtQ);

        mtmp.Reallocate(Nprevious_projection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, U_check,
               SeldonNoTrans, Q_sqrt_check, T(0), mtmp);
        mtmp_1 = FtQ;
        MltAdd(T(-1), FtU, mtmp, T(1), mtmp_1);
        mtmp.Reallocate(Nprojection_, Nmode_Q_);
        MltAdd(T(1), mtmp_1, Vinv, T(0), mtmp);
        MltAdd(T(-1), SeldonNoTrans, mtmp, SeldonTrans, mtmp_1, T(1), G_);

        // Adds $H_{t+1}^T R_{t+1}^{-1} H_{t+1}$ to 'G_'.
        Matrix<T, General, RowMajor> Ht_Rinv(Nprojection_, Nobservation_);
        MltAdd(T(1), SeldonTrans, H, SeldonNoTrans, R_inv_, T(0), Ht_Rinv);
        MltAdd(T(1), Ht_Rinv, H, T(1), G_);

        /*** Computes minimax state estimation ***/

        // Computes $z_t$.
        Add(T(1), e_, MFBz);
        for (int j = 0; j < Nstate_; j++)
            MFBz(j) *= D_tilde_inv_(j);

        z_.Reallocate(Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, MFBz, T(0), z_);

        vtmp.Reallocate(Nmode_Q_);
        MltAdd(T(1), SeldonTrans, Q_sqrt_check, MFBz, T(0), vtmp);
        MFBz.Clear();
        vtmp_1.Reallocate(Nmode_Q_);
        Mlt(IQtQinv, vtmp, vtmp_1);
        MltAdd(T(-1), SeldonNoTrans, FtQ, vtmp_1, T(1), z_);

        // Adds $H^T_{t+1} R_{t+1}^{-1} (y_{t+1} - \eta_{t+1})$.
        MltAdd(T(1), Ht_Rinv, y, T(1), z_);

        // Finally computes the minimax state, with $G_{t+1}^{-1} z_{t+1}$.
        state_ = z_;
        mtmp = G_;
        GetAndSolveLU(mtmp, state_);

        vtmp.Reallocate(Nstate_);
        MltAdd(T(1), SeldonTrans, projection_, state_, T(0), vtmp);
        model_.SetState(vtmp);

        if (inner_iteration_ == 0)
        {
            previous_projection_ = projection_;
            Nprevious_projection_ = previous_projection_.GetM();
        }

        inner_iteration_++;
        if (reduction_method_ != "none" && inner_iteration_ == Nsnapshot_ - 1)
        {
            inner_iteration_ = 0;
            mode_ = 0;
        }

        output_saver_.Save(G_, model_.GetTime(), "minimax_gain");
        output_saver_.Save(projection_, "projection");

        MessageHandler::Send(*this, "all", "::Propagation end");
    }


    //! Computes the full matrix of the tangent linear model.
    /*!
      \param[out] M the tangent linear model
    */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>
    ::ComputeTangentLinearModel(Matrix<T, General, RowMajor>& M)
    {
	Matrix<T, General, RowMajor> I(Nprevious_projection_,
                                       Nprevious_projection_);
	Vector<T> pi(Nstate_);
	Vector<T> Ei(Nprevious_projection_);

	for (int i = 0; i < Nprevious_projection_; i++)
	{
            Ei.Zero();
            Ei(i) = T(1);
	    // Takes the i-th vector from the basis, forming the projected
	    // subspace these vectors are stored as rows of the matrix
	    // 'projection_'.
	    MltAdd(T(1), SeldonTrans, previous_projection_, Ei, T(0), pi);
	    // Applies tangent linear model to the i-th vector from the basis.
            model_.ApplyTangentLinearOperator(pi);
	    SetCol(pi, i, M);
	}
    }


    //! Performs a step forward with snapshot recording.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::SnapshotRecording()
    {

        /*** Model ***/

        model_.Forward();
        model_state state;
        model_.GetState(state);
        SetCol(state, inner_iteration_, snapshot_);

        inner_iteration_++;

        /*** POD ***/

        if (inner_iteration_ == Nsnapshot_)
        {
            output_saver_.Save(snapshot_, "snapshot");

            GetSVD(snapshot_, singular_value_,
                   left_singular_vector_, right_singular_vector_);

            output_saver_.Save(singular_value_, "singular_value");
            output_saver_.Save(left_singular_vector_, "left_singular_vector");
            output_saver_.Save(right_singular_vector_,
                               "right_singular_vector");

            // Projection matrix.
            T total_energy = Norm1(singular_value_);
            T energy = 0;
            // Collects the main modes for the projection.
            int i;
            for (i = 0; i < Nprojection_max_
                     && energy <= (1. - acceptable_error_) * total_energy;
                 i++)
                energy += singular_value_(i);
            Nprojection_ = i;
            // Saving the previous projection because it is still useful for
            // the first step of reduced minimax.
            if (!first_sequence_)
            {
                previous_projection_ = projection_;
                Nprevious_projection_ = previous_projection_.GetM();
            }
            projection_.Reallocate(Nprojection_, Nstate_);
            for (i = 0; i < Nprojection_; i++)
                for (int j = 0; j < Nstate_; j++)
                    projection_(i, j) = left_singular_vector_(j, i);
            if (first_sequence_)
            {
                previous_projection_ = projection_;
                Nprevious_projection_ = previous_projection_.GetM();
            }

            Logger::Log(GetName(),
                        "Total energy and percentage retained with "
                        + to_str(Nprojection_) + " modes: "
                        + to_str(energy) + ", "
                        + to_str(energy / total_energy));

            inner_iteration_ = 0;
            mode_ = 1;
            iteration_ -= Nsnapshot_ - 1;
        }
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedMinimax<T, Model, ObservationManager>::HasFinished() const
    {
        // The condition 'mode_ == 0 && inner_iteration_ != Nsnapshot_' means
        // that a sequence of error computation has just ended. The simulation
        // should therefore end with a complete sequence of error computation.
        return (reduction_method_ == "none"
                || (mode_ == 0 && inner_iteration_ == 0))
            && model_.HasFinished();
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current projection matrix.
    /*!
      \return The current projection matrix.
    */
    template <class T, class Model, class ObservationManager>
    Matrix<T, General, RowMajor>&
    ReducedMinimax<T, Model, ObservationManager>::GetProjection()
    {
        return projection_;
    }


    //! Returns the previous projection matrix.
    /*!
      \return The previous projection matrix.
    */
    template <class T, class Model, class ObservationManager>
    Matrix<T, General, RowMajor>&
    ReducedMinimax<T, Model, ObservationManager>::GetPreviousProjection()
    {
        return previous_projection_;
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    const Model&
    ReducedMinimax<T, Model, ObservationManager>::GetModel() const
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class T, class Model, class ObservationManager>
    const ObservationManager&
    ReducedMinimax<T, Model, ObservationManager>
    ::GetObservationManager() const
    {
        return observation_manager_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string ReducedMinimax<T, Model, ObservationManager>::GetName() const
    {
        return "ReducedMinimax";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Message(string message)
    {
        model_state state;
        if (message.find("full_analysis") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, model_.GetTime(), "full_analysis");
        }
        else if (message.find("reduced_analysis") != string::npos)
            output_saver_.Save(state_, model_.GetTime(), "reduced_analysis");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDMINIMAX_CXX
#endif