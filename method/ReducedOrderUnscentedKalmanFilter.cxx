// Copyright (C) 2008-2010 INRIA
// Author(s): Marc Fragu, Philippe Moireau, Vivien Mallet
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_CXX

#include "ReducedOrderUnscentedKalmanFilter.hxx"

#include "seldon/vector/VectorCollection.cxx"

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
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::ReducedOrderUnscentedKalmanFilter(string configuration_file):
        configuration_file_(configuration_file)
    {

        /*** Initializations ***/

#if defined(VERDANDI_WITH_MPI)
        rank_ = MPI::COMM_WORLD.Get_rank();
        Nprocess_ = MPI::COMM_WORLD.Get_size();
        if (rank_ == 0)
        {
#endif
            MessageHandler::AddRecipient("model", model_,
                                         Model::StaticMessage);
            MessageHandler::AddRecipient("observation_manager",
                                         observation_manager_,
                                         ObservationManager::StaticMessage);
            MessageHandler::AddRecipient("driver", *this,
                                         ReducedOrderUnscentedKalmanFilter
                                         ::StaticMessage);
#if defined(VERDANDI_WITH_MPI)
        }
#endif

        /***************************
         * Reads the configuration *
         ***************************/


        Ops configuration(configuration_file_);
        configuration.SetPrefix("reduced_order_unscented_kalman_filter.");

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
        configuration.Set("data_assimilation.with_resampling",
                          with_resampling_);
        configuration.Set("data_assimilation.observation_error_variance",
                          "ops_in(v, {'matrix', 'matrix_inverse', 'vector',"
                          "'vector_inverse'})", observation_error_variance_);

        /*** Sigma-points ***/

        configuration.Set("sigma_point.type",
                          "ops_in(v, {'canonical', 'star', 'simplex'})",
                          sigma_point_type_);

#if defined(VERDANDI_WITH_MPI)

        /*** MPI ***/

        configuration.Set("mpi.algorithm", "ops_in(v, {0, 1, 2})",
                          algorithm_);
        configuration.Set("mpi.master_process_contribution",
                          master_process_contribution_);
        if (master_process_contribution_ < 0. ||
            master_process_contribution_ > 1.)
            throw "Contribution of process 0 should be in [0, 1] "
                "], but " + to_str(master_process_contribution_) +
                " was provided.";

#endif

        /*** Ouput saver ***/

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            configuration.
                SetPrefix("reduced_order_unscented_kalman_filter"
                          ".output_saver.");
            output_saver_.Initialize(configuration);
            output_saver_.Empty("state_forecast");
            output_saver_.Empty("state_analysis");

            /*** Logger and read configuration ***/

            configuration.SetPrefix("reduced_order_unscented_kalman_filter.");

            if (configuration.Exists("output.log"))
                Logger::SetFileName(configuration.Get<string>("output.log"));

            if (configuration.Exists("output.configuration"))
            {
                string output_configuration;
                configuration.Set("output.configuration",
                                  output_configuration);
                configuration.WriteLuaDefinition(output_configuration);
            }
#if defined(VERDANDI_WITH_MPI)
        }
        Logger::SetFileName("Verdandi_Proc" + to_str(rank_));
#endif

    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::~ReducedOrderUnscentedKalmanFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Initialize()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Initialize begin");

        /*** Initializations ***/

        model_.Initialize(model_configuration_file_);
        observation_manager_.Initialize(model_,
                                        observation_configuration_file_);
        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

#ifdef VERDANDI_DENSE
        model_.GetStateErrorVarianceSqrt(L_, U_);
        Copy(U_, U_inv_);
        GetInverse(U_inv_);
#else
        model_state_error_variance L_sparse, U_sparse;
        model_.GetStateErrorVarianceSqrt(L_sparse, U_sparse);
        Matrix<T, General, ArrayRowSparse> L_array, U_array;
        ConvertRowSparseToArrayRowSparse(L_sparse, L_array);
        ConvertRowSparseToArrayRowSparse(U_sparse, U_array);
        ConvertArrayRowSparseToDense(L_array, L_);
        ConvertArrayRowSparseToDense(U_array, U_);
        Copy(U_, U_inv_);
        GetInverse(U_inv_);
#endif

        Nreduced_ = U_.GetN();

        /*** Sigma-points ***/

        sigma_point_matrix V_trans;
        if (sigma_point_type_ == "canonical")
            ComputeCanonicalSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                       alpha_constant_);
        else if (sigma_point_type_ == "star")
            ComputeStarSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                  alpha_constant_);
        else if (sigma_point_type_ == "simplex")
            ComputeSimplexSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                     alpha_constant_);
        if (alpha_constant_)
            alpha_ = D_alpha_(0);

        Nsigma_point_ = V_trans.GetM();

        // Initializes transpose of I.
        sigma_point_matrix P_alpha_v(Nreduced_, Nreduced_);
        I_trans_.Reallocate(Nsigma_point_, Nreduced_);

        if (alpha_constant_)
        {
            MltAdd(T(alpha_), SeldonTrans, V_trans, SeldonNoTrans, V_trans,
                   T(0), P_alpha_v);
            GetInverse(P_alpha_v);
            GetCholesky(P_alpha_v);
            MltAdd(T(1), SeldonNoTrans, V_trans, SeldonTrans, P_alpha_v,
                   T(0), I_trans_);
        }
        else
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Initialize()", "Calculation not "
                                 "implemented for no constant alpha_i.");

        // Initializes D_v.
        D_v_.Reallocate(Nsigma_point_, Nsigma_point_);
        if (alpha_constant_)
            MltAdd(T(alpha_ * alpha_), SeldonNoTrans, I_trans_, SeldonTrans,
                   I_trans_, T(0), D_v_);
        else
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Initialize()", "Calculation not "
                                 "implemented for no constant alpha_i.");

#if defined(VERDANDI_WITH_MPI)

        /*** Local sigma-points ***/

        Nlocal_sigma_point_ = int(Nsigma_point_ / Nprocess_);
        int r = Nsigma_point_ % Nprocess_;
        if (Nprocess_ - rank_ - 1 < r)
            Nlocal_sigma_point_++;

        Nlocal_sigma_point_sum_.Reallocate(Nprocess_ + 1);
        Nlocal_sigma_point_sum_(0) = 0;
        for (int i = 0; i < Nprocess_; i++)
        {
            Nlocal_sigma_point_sum_(i + 1) = int(Nsigma_point_ / Nprocess_);
            if (Nprocess_ - i - 1 < r)
                Nlocal_sigma_point_sum_(i + 1) += 1;
            Nlocal_sigma_point_sum_(i + 1) += Nlocal_sigma_point_sum_(i);
        }

        int Nsigma_point_0;
        Nsigma_point_0 =  int(master_process_contribution_ *
                              Nlocal_sigma_point_sum_(1));
        int q;
        q = int((Nlocal_sigma_point_sum_(1) - Nsigma_point_0) /
                (Nprocess_ - 1));
        r = (Nlocal_sigma_point_sum_(1) - Nsigma_point_0) % (Nprocess_ - 1);
        for (int i = 1; i < Nprocess_; i++)
        {
            Nlocal_sigma_point_sum_(i + 1) -= Nlocal_sigma_point_sum_(1) -
                Nsigma_point_0;
            Nlocal_sigma_point_sum_(i + 1) += q * i;
            if (i <= r)
                Nlocal_sigma_point_sum_(i + 1) += i;
            else
                Nlocal_sigma_point_sum_(i + 1) += r;

        }
        Nlocal_sigma_point_sum_(1) = Nsigma_point_0;
        for (int i = Nlocal_sigma_point_sum_(rank_);
             i < Nlocal_sigma_point_sum_(rank_ + 1); i++)
            local_sigma_point_.PushBack(i);
        Nlocal_sigma_point_ = local_sigma_point_.GetM();

        I_trans_global_.Copy(I_trans_);
        I_trans_.Reallocate(Nlocal_sigma_point_, Nreduced_);
        sigma_point I_col;
        for (int i = 0; i < Nlocal_sigma_point_; i++)
        {
            GetRow(I_trans_global_, local_sigma_point_(i), I_col);
            SetRow(I_col, i, I_trans_);
        }

        /*** Local columns of covariance matrix ***/

        Nlocal_filtered_column_ = int(Nreduced_ / Nprocess_);
        r = Nreduced_ % Nprocess_;
        if (Nprocess_ - rank_ - 1 < r)
            Nlocal_filtered_column_++;

        Nlocal_filtered_column_sum_.Reallocate(Nprocess_ + 1);
        Nlocal_filtered_column_sum_(0) = 0;
        for (int i = 0; i < Nprocess_; i++)
        {
            Nlocal_filtered_column_sum_(i + 1) = int(Nreduced_ / Nprocess_);
            if (Nprocess_ - i - 1 < r)
                Nlocal_filtered_column_sum_(i + 1) += 1;
            Nlocal_filtered_column_sum_(i + 1)
                += Nlocal_filtered_column_sum_(i);
        }

        for (int i = Nlocal_filtered_column_sum_(rank_);
             i < Nlocal_filtered_column_sum_(rank_ + 1); i++)
            local_filtered_column_.PushBack(i);

#endif

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            MessageHandler::Send(*this, "model", "initial condition");
            MessageHandler::Send(*this, "driver", "initial condition");

            MessageHandler::Send(*this, "all", "::Initialize end");
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Initializes a step for the unscented Kalman filter.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::InitializeStep()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::InitializeStep begin");

        model_.InitializeStep();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Forward()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Forward begin");

        if (sigma_point_type_ == "simplex")
        {
#if defined(VERDANDI_WITH_MPI)

            if (algorithm_ == 0)
            {
                model_state x;
                model_.GetState(x);

                sigma_point_matrix X_i_trans_global;
                sigma_point x_col;
                X_i_trans_.Reallocate(Nlocal_sigma_point_, Nstate_);
                if (rank_ == 0)
                {
                    sigma_point_matrix tmp;
                    GetCholesky(U_inv_);
                    Copy(L_, tmp);
                    MltAdd(T(1), tmp, U_inv_, T(0), L_);

                    // Computes X_n^{(i)+}.
                    X_i_trans_global.Reallocate(Nsigma_point_, Nstate_);
                    for (int i = 0; i < Nsigma_point_; i++)
                        SetRow(x, i, X_i_trans_global);
                    MltAdd(T(1), SeldonNoTrans, I_trans_global_, SeldonTrans,
                           L_, T(1), X_i_trans_global);
                }

                int displacement[Nprocess_],  sendcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    sendcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nstate_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nstate_;
                }

                MPI::COMM_WORLD.Scatterv(X_i_trans_global.GetData(),
                                         sendcount, displacement, MPI::DOUBLE,
                                         X_i_trans_.GetData(),
                                         Nlocal_sigma_point_ * Nstate_,
                                         MPI::DOUBLE, 0);

                /*** Prediction ***/

                // Computes X_{n + 1}^-.
                x.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRow(X_i_trans_, i, x_col);
                    model_.ApplyOperator(x_col, i + 1 == Nlocal_sigma_point_,
                                         true);
                    Add(T(alpha_), x_col, x);
                    if (rank_ == 0)
                        SetRow(x_col, i, X_i_trans_global);
                    SetRow(x_col, i, X_i_trans_);
                }

                model_state working_vector(model_.GetNstate());
                working_vector.Fill(T(0));
                MPI::COMM_WORLD.Reduce(x.GetData(), working_vector.GetData(),
                                       x.GetM(), MPI::DOUBLE, MPI::SUM, 0);


                if (rank_ == 0)
                    model_.SetState(working_vector);

                MPI::COMM_WORLD.
                    Gatherv(X_i_trans_.GetData(),
                            Nlocal_sigma_point_ * Nstate_, MPI::DOUBLE,
                            X_i_trans_global.GetData(),  sendcount,
                            displacement, MPI::DOUBLE, 0);

                if (rank_ == 0)
                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_global,
                           SeldonNoTrans, I_trans_global_, T(0), L_);
            }
            else if (algorithm_ == 1 || algorithm_ == 2)
            {
                model_state x;
                model_.GetState(x);

                /*** Sampling ***/

                sigma_point_matrix tmp;
                GetCholesky(U_inv_);
                Copy(L_, tmp);
                MltAdd(T(1), tmp, U_inv_, T(0), L_);

                // Computes X_n^{(i)+}.
                X_i_trans_.Reallocate(Nlocal_sigma_point_, Nstate_);
                sigma_point x_col;
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                    SetRow(x, i, X_i_trans_);
                MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans, L_, T(1),
                       X_i_trans_);

                /*** Prediction ***/

                // Computes X_{n + 1}^-.
                x.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRow(X_i_trans_, i, x_col);
                    model_.ApplyOperator(x_col, i + 1 == Nlocal_sigma_point_,
                                         true);
                    Add(T(alpha_), x_col, x);
                    SetRow(x_col, i, X_i_trans_);
                }

                model_state working_vector(model_.GetNstate());
                working_vector.Fill(T(0));
                MPI::COMM_WORLD.Allreduce(x.GetData(),
                                          working_vector.GetData(), x.GetM(),
                                          MPI::DOUBLE, MPI::SUM);

                int displacement[Nprocess_],  recvcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    recvcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nstate_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nstate_;
                }

                // Computes L_{n + 1}.
                if (algorithm_ == 1)
                {
                    sigma_point_matrix
                        X_i_trans_global(Nsigma_point_, Nstate_);

                    MPI::COMM_WORLD.
                        Allgatherv(X_i_trans_.GetData(),
                                   Nlocal_sigma_point_ * Nstate_, MPI::DOUBLE,
                                   X_i_trans_global.GetData(),  recvcount,
                                   displacement, MPI::DOUBLE);

                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_global,
                           SeldonNoTrans, I_trans_global_, T(0), L_);
                }
                else
                {
                    sigma_point_matrix L_local(Nstate_, Nreduced_);
                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_, SeldonNoTrans,
                           I_trans_, T(0), L_local);

                    L_.Fill(T(0));
                    MPI::COMM_WORLD.
                        Allreduce(L_local.GetData(), L_.GetData(),
                                  L_.GetSize(), MPI::DOUBLE, MPI::SUM);
                }

                model_.SetState(working_vector);
            }
#else
            model_state x;
            model_.GetState(x);

            /*** Sampling ***/

            sigma_point_matrix tmp;
            GetCholesky(U_inv_);

            Copy(L_, tmp);
            MltAdd(T(1), tmp, U_inv_, T(0), L_);

            // Computes X_n^{(i)+}.
            X_i_trans_.Reallocate(Nsigma_point_, Nstate_);
            sigma_point x_col;
            for (int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, X_i_trans_);

            MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans, L_, T(1),
                   X_i_trans_);

            /*** Prediction ***/

            // Computes X_{n + 1}^-.
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRow(X_i_trans_, i, x_col);
                model_.ApplyOperator(x_col, i + 1 == Nsigma_point_, true);
                Add(T(alpha_), x_col, x);
                SetRow(x_col, i, X_i_trans_);
            }

            // Computes L_{n + 1}.
            MltAdd(T(alpha_), SeldonTrans, X_i_trans_, SeldonNoTrans,
                   I_trans_, T(0), L_);
            model_.SetState(x);

            /*** Resampling ***/

            if (with_resampling_)
            {
                for(int i = 0; i < Nsigma_point_; i++)
                    SetRow(x, i, X_i_trans_);
                MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans,
                       L_, T(1), X_i_trans_);
            }
#endif
        }
        else
        {
#if defined(VERDANDI_WITH_MPI)
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Forward()", "Parallel algorithm not "
                                 "implemented yet for the 'no"
                                 " simplex' cases.");
#else
            model_state x;
            model_.GetState(x);

            /*** Sampling ***/

            sigma_point_matrix tmp;
            GetCholesky(U_inv_);
            Copy(L_, tmp);
            MltAdd(T(1), tmp, U_inv_, T(0), L_);

            // Computes X_n^{(i)+}.
            X_i_trans_.Reallocate(Nsigma_point_, Nstate_);
            sigma_point x_col;
            for (int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, X_i_trans_);

            MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans, L_, T(1),
                   X_i_trans_);

            /*** Prediction ***/

            // Computes X_{n + 1}^-.
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRow(X_i_trans_, i, x_col);
                model_.ApplyOperator(x_col, i + 1 == Nsigma_point_, true);
                Add(T(alpha_), x_col, x);
                SetRow(x_col, i, X_i_trans_);
            }

            /*** Resampling with SVD ***/

            sigma_point_matrix M_trans(Nsigma_point_, Nstate_);
            for (int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, M_trans);
            Mlt(T(-1), M_trans);
            Add(T(1), X_i_trans_, M_trans);

            if (alpha_constant_)
                Mlt(sqrt(alpha_), M_trans);
            else
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Forward()", "Calculation not "
                                     "implemented for no constant alpha_i.");

            sigma_point_matrix G(Nsigma_point_, Nsigma_point_);
            MltAdd(T(1), SeldonNoTrans, M_trans, SeldonTrans, M_trans,
                   T(0), G);

            Vector<T> lambda;
            Matrix<T> U, V;
            GetSVD(G, lambda, U, V);
            U.Resize(Nsigma_point_, Nreduced_);

            sigma_point_matrix
                working_matrix_rr(Nsigma_point_, Nsigma_point_),
                working_matrix_rN(X_i_trans_);

            MltAdd(T(sqrt(alpha_)), SeldonNoTrans, U, SeldonTrans, I_trans_,
                   T(0), working_matrix_rr);

            for(int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, X_i_trans_);
            Add(T(-1), X_i_trans_, working_matrix_rN);

            MltAdd(T(1), SeldonTrans, working_matrix_rr, SeldonNoTrans,
                   working_matrix_rN, T(1), X_i_trans_);

            // Computes L_{n + 1}.
            MltAdd(T(alpha_), SeldonTrans, X_i_trans_, SeldonNoTrans,
                   I_trans_, T(0), L_);
            model_.SetState(x);
#endif
        }

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            MessageHandler::Send(*this, "model", "forecast");
            MessageHandler::Send(*this, "observation_manager", "forecast");
            MessageHandler::Send(*this, "driver", "forecast");

            MessageHandler::Send(*this, "all", "::Forward end");
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Analyze()
    {

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (!observation_manager_.HasObservation())
        {
#if defined(VERDANDI_WITH_MPI)
            if (rank_ == 0)
#endif
                MessageHandler::Send(*this, "all", "::Analyze end");
            return;
        }

        Nobservation_  = observation_manager_.GetNobservation();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            if (option_display_["show_time"])
                cout << "Performing Reduced Order UKF at time step ["
                     << model_.GetTime() << "]..." << endl;
#if defined(VERDANDI_WITH_MPI)
        }
#endif

        if (sigma_point_type_ == "simplex")
        {
#if defined(VERDANDI_WITH_MPI)
            if (algorithm_ == 0)
            {
                // Computes [HX_{n+1}^{*}].
                sigma_point_matrix Z_i_trans,
                    Z_i_trans_local(Nlocal_sigma_point_, Nobservation_);
                sigma_point x_col;
                observation z_col(Nobservation_), z_local(Nobservation_), z;
                z_local.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    SetRow(z_col, i, Z_i_trans_local);
                    Add(T(alpha_), z_col, z_local);
                    x_col.Nullify();
                }

                if (rank_ == 0)
                {
                    z.Reallocate(Nobservation_);
                    z.Fill(T(0));
                    Z_i_trans.Reallocate(Nsigma_point_, Nobservation_);
                }

                MPI::COMM_WORLD.
                    Reduce(z_local.GetData(), z.GetData(),
                           z_local.GetM(), MPI::DOUBLE, MPI::SUM, 0);

                int displacement[Nprocess_],  sendcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    sendcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nobservation_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nobservation_;
                }

                MPI::COMM_WORLD.
                    Gatherv(Z_i_trans_local.GetData(),
                            Nlocal_sigma_point_ * Nobservation_, MPI::DOUBLE,
                            Z_i_trans.GetData(),  sendcount, displacement,
                            MPI::DOUBLE, 0);

                model_state x;
                model_.GetState(x);
                if (rank_ == 0)
                {
                    sigma_point_matrix HL(Nobservation_, Nreduced_);
                    MltAdd(T(alpha_), SeldonTrans, Z_i_trans, SeldonNoTrans,
                           I_trans_global_, T(0), HL);
                    observation_error_variance R_inv;
                    sigma_point_matrix
                        working_matrix_po(Nreduced_, Nobservation_), tmp;

                    if (observation_error_variance_ == "matrix_inverse")
                        Copy(observation_manager_.GetErrorVarianceInverse(),
                             R_inv);
                    else if (observation_error_variance_ == "matrix")
                    {
                        Copy(observation_manager_.GetErrorVariance(),
                             R_inv);
                        GetInverse(R_inv);
                    }
                    else
                        throw ErrorUndefined(
                            "ReducedOrderUnscentedKalmanFilter::Analyse()",
                            "The parameter observation_error_variance = "
                            + observation_error_variance_ +
                            " is not implemented yet. Only 'matrix' and "
                            "'matrix_inverse' are implemented.");
#ifdef VERDANDI_DENSE
                    MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv,
                           T(0), working_matrix_po);
#else
                    sigma_point_matrix R_inv_dense;
                    ConvertRowSparseToDense(R_inv, R_inv_dense);
                    MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv_dense,
                           T(0), working_matrix_po);
#endif
                    U_inv_.SetIdentity();
                    MltAdd(T(1), working_matrix_po, HL, T(1), U_inv_);
                    GetInverse(U_inv_);

                    // Computes K;
                    sigma_point_matrix K(Nstate_, Nobservation_);
                    Copy(working_matrix_po, tmp);
                    MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
                    MltAdd(T(1), L_, tmp, T(0), K);

                    // Computes innovation.
                    observation innovation;
                    observation_manager_.GetObservation(innovation);
                    Add(T(-1), z, innovation);

                    // Updates.
                    MltAdd(T(1), K, innovation, T(1), x);
                    model_.SetState(x);
                }
            }
            else if (algorithm_ == 1 || algorithm_ == 2)
            {
                observation innovation(Nobservation_);
                MPI::Request send_request[Nprocess_ - 1], recv_request;
                if (rank_ == 0)
                {
                    observation_manager_.GetObservation(innovation);
                    for (int i = 1; i < Nprocess_; i++)
                        send_request[i] =
                            MPI::COMM_WORLD.
                            Isend(innovation.GetData(),
                                  Nobservation_, MPI::DOUBLE, i, 0);
                }
                else
                    recv_request =
                        MPI::COMM_WORLD.Irecv(innovation.GetData(),
                                              Nobservation_, MPI::DOUBLE, 0,
                                              MPI::ANY_TAG);
                model_state x;
                model_.GetState(x);
                // Computes [HX_{n+1}^{*}].
                sigma_point_matrix
                    Z_i_trans(Nlocal_sigma_point_, Nobservation_);
                sigma_point x_col;
                observation z_col, z(Nobservation_), z_local(Nobservation_);
                z_local.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    GetRowPointer(Z_i_trans, i, z_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    Add(T(alpha_), z_col, z_local);
                    x_col.Nullify();
                    z_col.Nullify();
                }

                sigma_point_matrix HL_local(Nobservation_, Nreduced_),
                    HL(Nobservation_, Nreduced_);
                MltAdd(T(alpha_), SeldonTrans, Z_i_trans, SeldonNoTrans,
                       I_trans_, T(0), HL_local);

                MPI::COMM_WORLD.Allreduce(HL_local.GetData(), HL.GetData(),
                                          HL.GetSize(), MPI::DOUBLE,
                                          MPI::SUM);

                observation_error_variance R_inv;
                sigma_point_matrix
                    working_matrix_po(Nreduced_, Nobservation_), tmp;

                if (observation_error_variance_ == "matrix_inverse")
                    Copy(observation_manager_.GetErrorVarianceInverse(),
                         R_inv);
                else if (observation_error_variance_ == "matrix")
                {
                    Copy(observation_manager_.GetErrorVariance(), R_inv);
                    GetInverse(R_inv);
                }
                else
                    throw ErrorUndefined(
                        "ReducedOrderUnscentedKalmanFilter::Analyse()",
                        "The parameter observation_error_variance = "
                        + observation_error_variance_ +
                        " is not implemented yet. Only 'matrix' and "
                        "'matrix_inverse' are implemented.");

                sigma_point_matrix K;
                if (algorithm_ == 1)
                {
#ifdef VERDANDI_DENSE
                    MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv,
                           T(0), working_matrix_po);
#else
                    sigma_point_matrix R_inv_dense;
                    ConvertRowSparseToDense(R_inv, R_inv_dense);
                    MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv_dense,
                           T(0), working_matrix_po);
#endif
                    U_inv_.SetIdentity();
                    MltAdd(T(1), working_matrix_po, HL, T(1), U_inv_);
                    GetInverse(U_inv_);

                    // Computes K;
                    K.Reallocate(Nstate_, Nobservation_);
                    Copy(working_matrix_po, tmp);
                    MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
                    MltAdd(T(1), L_, tmp, T(0), K);
                }
                else
                {
                    sigma_point_matrix working_matrix_po_local(
                        Nlocal_filtered_column_, Nobservation_);
                    sigma_point_matrix
                        U_local(Nlocal_filtered_column_, Nreduced_);

                    HL_local.
                        Reallocate(Nobservation_, Nlocal_filtered_column_);
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                    {
                        GetCol(HL, local_filtered_column_(i), z_col);
                        SetCol(z_col, i, HL_local);
                    }
#ifdef VERDANDI_DENSE
                    MltAdd(T(1), SeldonTrans, HL_local, SeldonNoTrans, R_inv,
                           T(0), working_matrix_po_local);
#else
                    sigma_point_matrix R_inv_dense;
                    ConvertRowSparseToDense(R_inv, R_inv_dense);
                    MltAdd(T(1), SeldonTrans, HL_local, SeldonNoTrans,
                           R_inv_dense, T(0), working_matrix_po_local);
#endif
                    U_local.Fill(T(0));
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                        U_local(i, local_filtered_column_(i)) = T(1);

                    MltAdd(T(1), working_matrix_po_local, HL, T(1), U_local);

                    int displacement[Nprocess_],  recvcount[Nprocess_];
                    for (int i = 0; i < Nprocess_; i++)
                    {
                        recvcount[i] = (Nlocal_filtered_column_sum_(i + 1)
                                        - Nlocal_filtered_column_sum_(i))
                            * Nreduced_;
                        displacement[i] = Nlocal_filtered_column_sum_(i)
                            * Nreduced_;
                    }

                    MPI::COMM_WORLD.
                        Allgatherv(U_local.GetData(),
                                   Nlocal_filtered_column_ * Nreduced_,
                                   MPI::DOUBLE, U_inv_.GetData(),  recvcount,
                                   displacement, MPI::DOUBLE);

                    for (int i = 0; i < Nprocess_; i++)
                    {
                        recvcount[i] = (Nlocal_filtered_column_sum_(i + 1)
                                        - Nlocal_filtered_column_sum_(i))
                            * Nobservation_;
                        displacement[i] = Nlocal_filtered_column_sum_(i)
                            * Nobservation_;
                    }

                    MPI::COMM_WORLD.
                        Allgatherv(working_matrix_po_local.GetData(),
                                   Nlocal_filtered_column_ * Nobservation_,
                                   MPI::DOUBLE, working_matrix_po.GetData(),
                                   recvcount, displacement, MPI::DOUBLE);

                    GetInverse(U_inv_);

                    // Computes K.
                    K.Reallocate(Nstate_, Nobservation_);
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                    {
                        GetRow(U_inv_, local_filtered_column_(i), z_col);
                        SetRow(z_col, i, U_local);
                    }

                    MltAdd(T(1), U_local, working_matrix_po, T(0),
                           working_matrix_po_local);

                    MPI::COMM_WORLD.
                        Allgatherv(working_matrix_po_local.GetData(),
                                   Nlocal_filtered_column_ * Nobservation_,
                                   MPI::DOUBLE, working_matrix_po.GetData(),
                                   recvcount, displacement, MPI::DOUBLE);

                    MltAdd(T(1), L_, working_matrix_po, T(0), K);
                }

                // Computes innovation.
                if (rank_ == 0)
                    MPI::Request::Waitall(Nprocess_ - 1, send_request);
                else
                    recv_request.Wait();

                MPI::COMM_WORLD.Allreduce(z_local.GetData(), z.GetData(),
                                          Nobservation_, MPI::DOUBLE,
                                          MPI::SUM);
                Add(T(-1), z, innovation);

                // Updates.
                MltAdd(T(1), K, innovation, T(1), x);
                model_.SetState(x);
            }
#else
            model_state x;
            model_.GetState(x);

            // Computes [HX_{n+1}^{*}].
            sigma_point_matrix Z_i_trans(Nsigma_point_, Nobservation_);
            sigma_point x_col;
            observation z_col, z(Nobservation_);
            z.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRowPointer(X_i_trans_, i, x_col);
                GetRowPointer(Z_i_trans, i, z_col);
                observation_manager_.ApplyOperator(x_col, z_col);
                Add(T(alpha_), z_col, z);
                x_col.Nullify();
                z_col.Nullify();
            }
            sigma_point_matrix HL(Nobservation_, Nreduced_);
            MltAdd(T(alpha_), SeldonTrans, Z_i_trans, SeldonNoTrans, I_trans_,
                   T(0), HL);

            observation_error_variance R_inv;
            sigma_point_matrix working_matrix_po(Nreduced_, Nobservation_),
                tmp;

            if (observation_error_variance_ == "matrix_inverse")
                Copy(observation_manager_.GetErrorVarianceInverse(), R_inv);
            else if (observation_error_variance_ == "matrix")
            {
                Copy(observation_manager_.GetErrorVariance(), R_inv);
                GetInverse(R_inv);
            }
            else
                throw ErrorUndefined(
                    "ReducedOrderUnscentedKalmanFilter::Analyse()",
                    "The parameter observation_error_variance = "
                    + observation_error_variance_ +
                    " is not implemented yet. Only 'matrix' and "
                    "'matrix_inverse' are implemented.");
#ifdef VERDANDI_DENSE
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv,
                   T(0), working_matrix_po);
#else
            sigma_point_matrix R_inv_dense;
            ConvertRowSparseToDense(R_inv, R_inv_dense);
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv_dense,
                   T(0), working_matrix_po);
#endif
            U_inv_.SetIdentity();
            MltAdd(T(1), working_matrix_po, HL, T(1), U_inv_);
            GetInverse(U_inv_);

            // Computes K;
            sigma_point_matrix K(Nstate_, Nobservation_);
            Copy(working_matrix_po, tmp);
            MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
            MltAdd(T(1), L_, tmp, T(0), K);

            // Computes innovation.
            observation innovation;
            observation_manager_.GetObservation(innovation);
            Add(T(-1), z, innovation);

            // Updates.
            MltAdd(T(1), K, innovation, T(1), x);
            model_.SetState(x);
#endif
        }
        else
        {
#if defined(VERDANDI_WITH_MPI)
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Analyse()", "Parallel algorithm not"
                                 " implemented yet for the 'no"
                                 " simplex' cases.");
#else
            model_state x;
            model_.GetState(x);

            // Computes [HX_{n+1}^{*}].
            sigma_point_matrix Z_i_trans(Nsigma_point_, Nobservation_);
            sigma_point x_col;
            observation z_col, z(Nobservation_);
            z.Fill(T(0));
            if (alpha_constant_)
            {
                for (int i = 0; i < Nsigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    GetRowPointer(Z_i_trans, i, z_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    Add(T(alpha_), z_col, z);
                    x_col.Nullify();
                    z_col.Nullify();
                }
            }
            else
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Analyse()", "Calculation not "
                                     "implemented for no constant alpha_i.");

            // Computes [Z] = [HX_{n+1}^{*} - E(HX_{n+1}^{*})].
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRowPointer(Z_i_trans, i, z_col);
                Add(T(-1), z, z_col);
                z_col.Nullify();
            }

            observation_error_variance R_inv;
            sigma_point_matrix
                working_matrix_ro(Nsigma_point_, Nobservation_),
                D_m(Nsigma_point_, Nsigma_point_);
            sigma_point_matrix HL;

            if (observation_error_variance_ == "matrix_inverse")
                Copy(observation_manager_.GetErrorVarianceInverse(), R_inv);
            else if (observation_error_variance_ == "matrix")
            {
                Copy(observation_manager_.GetErrorVariance(), R_inv);
                GetInverse(R_inv);
            }
            else
                throw ErrorUndefined(
                    "ReducedOrderUnscentedKalmanFilter::Analyse()",
                    "The parameter observation_error_variance = "
                    + observation_error_variance_ +
                    " is not implemented yet. Only 'matrix' and "
                    "'matrix_inverse' are implemented.");

            // Computes D_m.
            MltAdd(T(1), Z_i_trans, R_inv, T(0), working_matrix_ro);
            MltAdd(T(1), SeldonNoTrans, working_matrix_ro, SeldonTrans,
                   Z_i_trans, T(0), D_m);

            // Computes U_{n+1}.
            sigma_point_matrix
                working_matrix_rp(Nsigma_point_, Nreduced_),
                working_matrix_rr(Nsigma_point_, Nsigma_point_),
                working_matrix_rr2(Nsigma_point_, Nsigma_point_),
                working_matrix_rr3(Nsigma_point_, Nsigma_point_);

            Copy(D_v_, working_matrix_rr);
            Mlt(T(-1), working_matrix_rr);
            if (alpha_constant_)
            {
                for(int i = 0; i < Nsigma_point_; i++ )
                    working_matrix_rr(i, i) += alpha_;
                MltAdd(T(1), D_m, working_matrix_rr, T(0),
                       working_matrix_rr2);
                for(int i = 0; i < Nsigma_point_; i++ )
                    working_matrix_rr2(i, i) += 1;
                GetInverse(working_matrix_rr2);
                MltAdd(T(1), working_matrix_rr2, D_m, T(0),
                       working_matrix_rr);
                MltAdd(T(alpha_), working_matrix_rr, I_trans_, T(0),
                       working_matrix_rp);
                U_.SetIdentity();
                MltAdd(T(alpha_), SeldonTrans, I_trans_, SeldonNoTrans,
                       working_matrix_rp, T(1), U_);

                Copy(U_, U_inv_);
                GetInverse(U_inv_);

                // Computes {HL}_{n+1}.
                HL.Reallocate(Nobservation_, Nreduced_);
                working_matrix_rr2.SetIdentity();
                MltAdd(T(1), D_v_, working_matrix_rr, T(1),
                       working_matrix_rr2);

                working_matrix_rr.SetIdentity();
                Add(T(alpha_), D_m, working_matrix_rr);
                GetInverse(working_matrix_rr);

                MltAdd(T(1), working_matrix_rr, working_matrix_rr2, T(0),
                       working_matrix_rr3);

                MltAdd(T(alpha_), working_matrix_rr3, I_trans_, T(0),
                       working_matrix_rp);

                MltAdd(T(1), SeldonTrans, Z_i_trans, SeldonNoTrans,
                       working_matrix_rp, T(0), HL);
            }
            else
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Analyse()", "Calculation not "
                                     "implemented for no constant alpha_i.");


            // Computes K;
            sigma_point_matrix K(Nstate_, Nobservation_),
                working_matrix_po(Nreduced_, Nobservation_),
                working_matrix_po2(Nreduced_, Nobservation_);

#ifdef VERDANDI_DENSE
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv, T(0),
                   working_matrix_po);
            MltAdd(T(1), U_inv_, working_matrix_po, T(0), working_matrix_po2);
            MltAdd(T(1), L_, working_matrix_po2, T(0), K);
#else
            sigma_point_matrix R_inv_dense;
            ConvertRowSparseToDense(R_inv, R_inv_dense);
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans, R_inv_dense, T(0),
                   working_matrix_po);
            MltAdd(T(1), U_inv_, working_matrix_po, T(0), working_matrix_po2);
            MltAdd(T(1), L_, working_matrix_po2, T(0), K);
#endif

            // Computes innovation.
            observation innovation;
            observation_manager_.GetObservation(innovation);
            Add(T(-1), z, innovation);

            // Updates.
            MltAdd(T(1), K, innovation, T(1), x);
            model_.SetState(x);

#endif
        }

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif

            if (option_display_["show_time"])
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");

            MessageHandler::Send(*this, "all", "::Analyze end");

#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::GetName() const
    {
        return "ReducedOrderUnscentedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Message(string message)
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
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
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_CXX
#endif
