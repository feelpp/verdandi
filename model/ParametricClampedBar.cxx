// Copyright (C) 2008-2009 INRIA
// Author(s): Dominique Chapelle, Philippe Moireau, Marc Fragu
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


#ifndef VERDANDI_FILE_MODEL_PARAMETRICCLAMPEDBAR_CXX


#include "ParametricClampedBar.hxx"

#include "OutputSaver.cxx"

#include "seldon/vector/VectorCollection.cxx"

namespace Verdandi
{


    ///////////////////
    // STATIC FIELDS //
    ///////////////////


    template <class T>
    const double ParametricClampedBar<T>::Pi_ = 3.141592653589793238462;


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    template <class T>
    ParametricClampedBar<T>::ParametricClampedBar(): time_(0.)
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    ParametricClampedBar<T>::ParametricClampedBar(string configuration_file):
        time_(0.)
    {
    }


    //! Destructor.
    template <class T>
    ParametricClampedBar<T>::~ParametricClampedBar()
    {
        x_.Nullify();
        x_full_.Nullify();
    }


    ////////////////
    // INITIALIZE //
    ////////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ParametricClampedBar<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        Ops configuration(configuration_file);

        configuration.SetPrefix("parametric_clamped_bar.domain.");
        configuration.Set("bar_length", bar_length_);
        configuration.Set("Nx", Nx_);
        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("final_time", final_time_);
        configuration.SetPrefix("parametric_clamped_bar.error_statistics.");
        configuration.Set("state_error_variance", "v >= 0",
                          state_error_variance_value_);
        configuration.Set("state_error_scale", "v > 0",
                          Balgovind_scale_background_);
        configuration.SetPrefix("parametric_clamped_bar.physics.");
        configuration.Set("Young_modulus", Young_modulus_);
        configuration.Set("mass_density", mass_density_);

        configuration.Set("theta_force", theta_force_);
        Ntheta_force_ = theta_force_.GetSize();
        Ndof_ = Nx_ + 1;
        BuildRegionIndex(Ndof_, Ntheta_force_, theta_force_index_);
        configuration.Set("theta_stiffness", theta_stiffness_);
        Ntheta_stiffness_ = theta_stiffness_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_stiffness_, theta_stiffness_index_);
        configuration.Set("theta_mass", theta_mass_);
        Ntheta_mass_ = theta_mass_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_mass_, theta_mass_index_);
        configuration.Set("theta_damp", theta_damp_);
        Ntheta_damp_ = theta_damp_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_damp_, theta_damp_index_);

        configuration.Set("alpha", alpha_);
        configuration.Set("beta", beta_);
        configuration.SetPrefix("parametric_clamped_bar.");
        vector<string> stable;
        configuration.Set("state","ops_in(v, {'displacement', 'velocity', "
                          "'theta_force', 'theta_stiffness', 'theta_mass', "
                          "'theta_damp'})",
                          stable);
        if (stable.size() <= 0)
            throw ErrorArgument("void ParametricClampedBar<T>::Initialize"
                                "(string configuration_file)", "The name of"
                                " the underlying state vector (state) are "
                                "not defined.");
        for (unsigned int i = 0; i < stable.size(); i++)
            stable_.insert(stable[i]);

        configuration.Set("reduced_state", reduced_);

        /*** Ouput saver ***/

        output_saver_.Initialize(configuration_file,
                                 "parametric_clamped_bar.output_saver.");
        output_saver_.Empty("disp_0");
        output_saver_.Empty("velo_0");

        /*** Allocation ***/

        Delta_x_ = bar_length_ / Nx_;
        AllocateSparseMatrix();

        /*** Elementary mass matrix construction ***/

        mass_FEM_matrix_.Reallocate(2, 2);
        T mass_lin = T(mass_density_ * Delta_x_ / 3);
        mass_FEM_matrix_(0, 0) = mass_lin;
        mass_FEM_matrix_(1, 1) = mass_lin;
        mass_FEM_matrix_(0, 1) = mass_lin / 2;

        /*** Elementary stifness matrix construction ***/

        stiffness_FEM_matrix_.Reallocate(2, 2);
        T stiff_lin = T(Young_modulus_ / Delta_x_);
        stiffness_FEM_matrix_(0, 0) = stiff_lin;
        stiffness_FEM_matrix_(1, 1) = stiff_lin;
        stiffness_FEM_matrix_(0, 1) = -stiff_lin;

        /*** State initialization ***/

        disp_0_.Reallocate(Ndof_);
        velo_0_.Reallocate(Ndof_);
        force_.Reallocate(Ndof_);
        disp_0_.Fill(T(0));
        velo_0_.Fill(T(0));
        force_.Fill(T(1));

        for (int i = 0; i < Ndof_; i++)
            disp_0_(i) = T(i) / T(Ndof_ - 1);

        // Initialize full vector collection.
        state working_vector;
        working_vector.SetData(Ndof_ - 1, disp_0_.GetData() + 1);
        x_full_.AddVector(working_vector, "displacement");
        working_vector.Nullify();
        working_vector.SetData(Ndof_ - 1, velo_0_.GetData() + 1);
        x_full_.AddVector(working_vector, "velocity");
        working_vector.Nullify();
        x_full_.AddVector(theta_force_, "theta_force");
        x_full_.AddVector(theta_stiffness_, "theta_stiffness");
        x_full_.AddVector(theta_mass_, "theta_mass");
        x_full_.AddVector(theta_damp_, "theta_damp");

        // Initialize vector collection.
        if (stable_.find("displacement") != stable_.end())
            x_.AddVector(x_full_.GetVector("displacement"), "displacement");
        if (stable_.find("velocity") != stable_.end())
            x_.AddVector(x_full_.GetVector("velocity"), "velocity");
        if (stable_.find("theta_force") != stable_.end())
            x_.AddVector(theta_force_, "theta_force");
        if (stable_.find("theta_stiffness") != stable_.end())
            x_.AddVector(theta_stiffness_, "theta_stiffness");
        if (stable_.find("theta_mass") != stable_.end())
            x_.AddVector(theta_mass_, "theta_mass");
        if (stable_.find("theta_damp") != stable_.end())
            x_.AddVector(theta_damp_, "theta_damp");

        Nstate_ = x_.GetM();

#ifdef VERDANDI_STATE_ERROR_SPARSE
        build_diagonal_sparse_matrix(GetNstate(),
                                     state_error_variance_value_,
                                     state_error_variance_);
#else
        state_error_variance_.Reallocate(GetNstate(), GetNstate());
        state_error_variance_.SetIdentity();
        Mlt(T(state_error_variance_value_), state_error_variance_);
#endif
    }


    //! Initializes the first time step for the model.
    template <class T>
    void ParametricClampedBar<T>::InitializeFirstStep()
    {
    }


    //! Initializes the current time step for the model.
    template <class T>
    void ParametricClampedBar<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    /*
      \param[in] update_force Boolean to indicate if the force has to be
      updated or not.
    */
    template <class T>
    void ParametricClampedBar<T>::Forward(bool update_force)
    {
        /*** Update time ***/

        time_ += Delta_t_;

        /*** Right hand side ***/

        if (update_force)
        {
            AssembleMassMatrix(theta_force_, theta_force_index_);
            state ones(Ndof_);
            ones.Fill(T(1));
            MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                   mass_matrix_, ones, T(0), force_);
            force_(0) = T(0);
        }

        /*** Assembling process ***/

        AssembleMassMatrix(theta_mass_, theta_mass_index_);

        AssembleDampMatrix();

        AssembleNewMarkMatrix0();
        AssembleNewMarkMatrix1();

#if defined(VERDANDI_WITH_DIRECT_SOLVER)
#if defined(SELDON_WITH_UMFPACK)
        MatrixUmfPack<T> mat_lu;
#elif defined(SELDON_WITH_SUPERLU)
        MatrixSuperLU<T> mat_lu;
#elif defined(SELDON_WITH_MUMPS)
        MatrixMumps<T> mat_lu;
#endif
        GetLU(Newmark_matrix_1_, mat_lu, true);
#endif

        // FEM Vector (disp 1).
        state disp_1(Ndof_);
        // FEM Vector (velo 1).
        state velo_1(Ndof_);
        velo_1.Fill(T(0));

        MltAdd(T(2) / Delta_t_, mass_matrix_, velo_0_, T(1), force_);
        MltAdd(T(1), Newmark_matrix_0_, disp_0_, T(1), force_);

        // Dirichlet conditions .
        force_(0) = T(0);

#ifdef VERDANDI_WITH_DIRECT_SOLVER
        Vector<T> tmp(Ndof_);
        tmp = force_;
        SolveLU(mat_lu, tmp);
        disp_1 = tmp;
#else
        // Initialization of the Gmres parameters.
        int nb_max_iter = 1000;
        double tolerance = 1e-6;
        Iteration<double> iter(nb_max_iter, tolerance);
        Preconditioner_Base precond;
        // No preconditioning.
        iter.SetRestart(5);
        iter.HideMessages();
        Gmres(Newmark_matrix_1_, disp_1, force_, precond, iter);
#endif

        velo_1.Fill(T(0));
        Add(T(-1), velo_0_, velo_1);
        Add(T(2) / Delta_t_, disp_1, velo_1);
        Add(T(-2) / Delta_t_, disp_0_, velo_1);

        // Update.
        disp_0_ = disp_1;
        velo_0_ = velo_1;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool ParametricClampedBar<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Saves the simulated data.
    /*! It saves the displacement 'disp_0_' and  the velocity 'velo_0_'.
     */
    template <class T>
    void ParametricClampedBar<T>::Save()
    {
        output_saver_.Save(disp_0_, time_, "disp_0");
        output_saver_.Save(velo_0_, time_, "velo_0");
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the model to a given state vector.
    /*!
      \param[in,out] x on entry, the state vector to which the model is
      applied; on exit, the state vector after the model is applied.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
    */
    template <class T>
    void ParametricClampedBar<T>::ApplyOperator(state& x, bool forward,
                                                bool preserve_state,
                                                bool update_force)
    {
        double saved_time = 0;
        state saved_state;
        if (!forward)
            saved_time = GetTime();

        if (preserve_state)
            GetFullState(saved_state);

        SetState(x);
        Forward(update_force);
        GetState(x);

        if (!forward)
            SetTime(saved_time);

        if (preserve_state)
            SetFullState(saved_state);
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in,out] increment the increment.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::ApplyTangentLinearOperator(state& increment)
    {
        double saved_time = 0;
        state saved_state;
        saved_time = GetTime();
        GetFullState(saved_state);

        // Saves state parameters.
        Vector<T> theta_force(Ntheta_force_), theta_damp(Ntheta_damp_),
            theta_stiffness(Ntheta_stiffness_), theta_mass(Ntheta_mass_);
        if (stable_.find("theta_force") != stable_.end())
            Copy(theta_force_, theta_force);
        if (stable_.find("theta_stiffness") != stable_.end())
            Copy(theta_stiffness_, theta_stiffness);
        if (stable_.find("theta_mass") != stable_.end())
            Copy(theta_mass_, theta_mass);
        if (stable_.find("theta_damp") != stable_.end())
            Copy(theta_damp_, theta_damp);

        // Computes x_{n+1} + x_n and v_{n+1} - v_n.
        state disp(Ndof_), velo(Ndof_);

        Copy(disp_0_, disp);
        Copy(velo_0_, velo);

        Forward(true);

        Add(T(1), disp_0_, disp);
        Mlt(T(-1), velo);
        Add(T(1), velo_0_, velo);

        /*** Right hand side ***/

        SetTime(saved_time);
        SetState(increment);

        // Saves increment parameters.
        Vector<T> delta_theta_force(Ntheta_force_),
            delta_theta_damp(Ntheta_damp_),
            delta_theta_stiffness(Ntheta_stiffness_),
            delta_theta_mass(Ntheta_mass_);
        if (stable_.find("theta_force") != stable_.end())
            Copy(theta_force_, delta_theta_force);
        if (stable_.find("theta_stiffness") != stable_.end())
            Copy(theta_stiffness_, delta_theta_stiffness);
        if (stable_.find("theta_mass") != stable_.end())
            Copy(theta_mass_, delta_theta_mass);
        if (stable_.find("theta_damp") != stable_.end())
            Copy(theta_damp_, delta_theta_damp);

        state force(Ndof_);
        force.Fill(0);
        // Assemble F.
        if (stable_.find("theta_force") != stable_.end())
        {
            AssembleMassMatrix(theta_force_, theta_force_index_);
            state ones(Ndof_);
            ones.Fill(T(1));
            MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                   mass_matrix_, ones, T(0), force);
            force(0) = T(0);
        }

        // Assemble K.
        if (stable_.find("theta_stiffness") != stable_.end())
        {
            Fill(T(0), Newmark_matrix_0_);
            for (int i = 0; i < Nx_; i++)
            {
                Newmark_matrix_0_.Val(i, i) += stiffness_FEM_matrix_(0, 0)
                    * theta_stiffness_(theta_stiffness_index_(i));
                Newmark_matrix_0_.Val(i + 1, i + 1) +=
                    stiffness_FEM_matrix_(1, 1)
                    * theta_stiffness_(theta_stiffness_index_(i));
                Newmark_matrix_0_.Val(i, i + 1) += stiffness_FEM_matrix_(0, 1)
                    * theta_stiffness_(theta_stiffness_index_(i));
            }
            MltAdd(T(-0.5), Newmark_matrix_0_, disp, T(1), force);
            force(0) = T(0);
        }

        // Assemble M.
        if (stable_.find("theta_mass") != stable_.end())
        {
            AssembleMassMatrix(theta_mass_, theta_mass_index_);
            MltAdd(T(T(-1) / Delta_t_), mass_matrix_, velo, T(1), force);
            force(0) = T(0);
        }

        Copy(force, force_);

        /*** Computes delta_y ***/

        // Sets state parameters.
        if (stable_.find("theta_force") != stable_.end())
            Copy(theta_force, theta_force_);
        if (stable_.find("theta_stiffness") != stable_.end())
            Copy(theta_stiffness, theta_stiffness_);
        if (stable_.find("theta_mass") != stable_.end())
            Copy(theta_mass, theta_mass_);
        if (stable_.find("theta_damp") != stable_.end())
            Copy(theta_damp, theta_damp_);

        Forward(false);

        // Sets increment parameters.
        if (stable_.find("theta_force") != stable_.end())
            Copy(delta_theta_force, theta_force_);
        if (stable_.find("theta_stiffness") != stable_.end())
            Copy(delta_theta_stiffness, theta_stiffness_);
        if (stable_.find("theta_mass") != stable_.end())
            Copy(delta_theta_mass, theta_mass_);
        if (stable_.find("theta_damp") != stable_.end())
            Copy(delta_theta_damp, theta_damp_);

        GetState(increment);

        SetTime(saved_time);
        SetFullState(saved_state);
    }


    //! Gets the matrix of the tangent linear model.
    /*!
      \param[out] A the matrix of the tangent linear model.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::GetTangentLinearOperator(tangent_linear_operator& A) const
    {
        throw ErrorUndefined("ParametricClampedBar::GetTangentLinearOperator"
                             "(tangent_linear_operator& A) const");
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T>
    double ParametricClampedBar<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the time of the model to a given time.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void ParametricClampedBar<T>::SetTime(double time)
    {
        time_= time;
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    template <class T>
    int ParametricClampedBar<T>::GetNstate() const
    {
        return Nstate_;
    }


    //! Provides the reduced state vector.
    /*!
      \param[out] state the reduced state vector.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::GetState(state& x) const
    {
        x.Reallocate(x_.GetM());
        for (int i = 0; i < x_.GetM(); i++)
            x(i) = x_(i);
    }


    //! Sets the reduced state vector.
    /*! Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::SetState(state& x)
    {
        if (x_.GetM() != x.GetM())
            throw ErrorProcessing("ParametricClampedBar::SetState()",
                                  "Operation not permitted:\n x_ is a vector"
                                  " of length " + to_str(x_.GetM()) +
                                  ";\n x is a vector of length "
                                  + to_str(x.GetM()) + ".");
        for (int i = 0; i < x_.GetM(); i++)
            x_(i) = x(i);
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    void ParametricClampedBar<T>::GetFullState(state& x) const
    {
        x.Reallocate(x_full_.GetM());
        for (int i = 0; i < x_full_.GetM(); i++)
            x(i) = x_full_(i);
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    template <class T>
    void ParametricClampedBar<T>::SetFullState(const state& x)
    {
        if (x_full_.GetM() != x.GetM())
            throw ErrorProcessing("ParametricClampedBar::SetState()",
                                  "Operation not permitted:\n x_full_ is a "
                                  "vector of length " + to_str(x_full_.GetM())
                                  + ";\n x is a vector of length "
                                  + to_str(x.GetM()) + ".");
        for (int i = 0; i < x_.GetM(); i++)
            x_full_(i) = x(i);
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] error_covariance_row the value of row number \a row.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::GetStateErrorVarianceRow(int row, state_error_variance_row&
                               state_error_variance_row)
    {
        GetRow(state_error_variance_, row, state_error_variance_row);
    }


    //! Returns the background error covariance matrix (\f$B\f$).
    /*! Returns the background error covariance matrix (\f$B\f$).
      \return The matrix of the background error covariance.
    */
    template <class T>
    typename ParametricClampedBar<T>::state_error_variance&
    ParametricClampedBar<T>::GetStateErrorVariance()
    {
        return state_error_variance_;
    }


    //! Returns the background error covariance matrix (\f$B\f$).
    /*! Returns the background error covariance matrix (\f$B\f$).
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename ParametricClampedBar<T>::state_error_variance&
    ParametricClampedBar<T>::GetStateErrorVariance() const
    {
        return state_error_variance_;
    }


    /*! Returns a decomposition of the state error covariance matrix (\f$B\f$)
      as a product \f$LUL^T\f$.
    */
    /*!
      \param[out] L the matrix \f$L\f$.
      \param[out] U the matrix \f$U\f$.
    */
    template <class T>
    void ParametricClampedBar<T>::GetStateErrorVarianceSqrt(
        state_error_variance& L,
        state_error_variance& U)
    {
        int Nreduced = 0;
        for (unsigned int i = 0; i < reduced_.size(); i++)
            Nreduced += x_.GetVector(reduced_[i]).GetSize();
#ifndef VERDANDI_STATE_ERROR_SPARSE
        // Initializes L.
        L.Reallocate(Nstate_, Nreduced);
        L.Fill(T(0));
        for (unsigned int i = 0, l = 0; i < reduced_.size(); i++)
            for(int k = x_.GetIndex(reduced_[i]);
                k < x_.GetIndex(reduced_[i]) +
                    x_.GetVector(reduced_[i]).GetSize(); k++)
                L(k, l++) = 1;

        // Initializes U.
        U.Reallocate(Nreduced,  Nreduced);
        U.Fill(T(0));
        for (int i = 0; i < Nreduced; i++)
            U(i, i) = T(T(1) / state_error_variance_value_);
#else
        // Initializes L.
        Matrix<T, General, ArrayRowSparse> L_array(Nstate_, Nreduced);
        for (unsigned int i = 0, l = 0; i < reduced_.size(); i++)
            for(int k = x_.GetIndex(reduced_[i]);
                k < x_.GetIndex(reduced_[i]) +
                    x_.GetVector(reduced_[i]).GetSize(); k++)
                L_array(k, l++) = 1;
        Copy(L_array, L);

        // Initializes U.
        Matrix<T, General, ArrayRowSparse> U_array(Nreduced,  Nreduced);
        for (int i = 0; i < Nreduced; i++)
            U_array(i, i) = T(T(1) / state_error_variance_value_);
        Copy(U_array, U);
#endif
    }


    //! Checks if the error covariance matrix is sparse.
    /*!
      \return True if there is a sparse error matrix, false otherwise.
    */
    template <class T>
    bool ParametricClampedBar<T>::IsErrorSparse() const
    {
#ifdef VERDANDI_STATE_ERROR_SPARSE
        return true;
#else
        return false;
#endif

    }


    //! Returns the name of the class.
    template <class T>
    string ParametricClampedBar<T>::GetName() const
    {
        return "ParametricClampedBar";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void ParametricClampedBar<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


    /*! Build a vector indexed by the points of the bar to indicate
      which region each point belongs to.
    */
    /*
      \param[in] N the size of the bar.
      \param[in] Nregion the number of region.
      \param[out] index_vector the vector of indexes.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::BuildRegionIndex(int N, int Nregion, Vector<int>& index_vector)
    {
        index_vector.Reallocate(N);
        for(int i = 0; i < N; i++)
            index_vector(i) = i % Nregion;
        Sort(index_vector);
    }


    //! Assembles the Mass matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void ParametricClampedBar<T>
    ::AssembleMassMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        Fill(T(0), mass_matrix_);
        for (int i = 0; i < Nx_; i++)
        {
            // Mass Matrix.
            mass_matrix_.Val(i, i) += theta(theta_index(i))
                * mass_FEM_matrix_(0, 0);
            mass_matrix_.Val(i + 1, i + 1) += theta(theta_index(i))
                * mass_FEM_matrix_(1, 1);
            mass_matrix_.Val(i, i + 1) += theta(theta_index(i))
                * mass_FEM_matrix_(0, 1);
        }

        // Boundary condition by pseudo-elimination.
        mass_matrix_.Val(0, 0) = 1;
        mass_matrix_.Val(0, 1) = 0;
        mass_matrix_.Val(1, 0) = 0;
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void ParametricClampedBar<T>
    ::AssembleNewMarkMatrix0()
    {
        Fill(T(0), Newmark_matrix_0_);
        for (int i = 0; i < Nx_; i++)
        {
            // Newmark's matrix at time n.
            Newmark_matrix_0_.Val(i, i) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(0, 0) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i, i) / Delta_t_ +
                -0.5 * stiffness_FEM_matrix_(0, 0)
                * theta_stiffness_(theta_stiffness_index_(i));
            Newmark_matrix_0_.Val(i + 1, i + 1) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(1, 1) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i + 1, i + 1) / Delta_t_ +
                -0.5 * stiffness_FEM_matrix_(1, 1)
                * theta_stiffness_(theta_stiffness_index_(i));
            Newmark_matrix_0_.Val(i, i + 1) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(0, 1) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i, i + 1) / Delta_t_ +
                -0.5 * stiffness_FEM_matrix_(0, 1)
                * theta_stiffness_(theta_stiffness_index_(i));
        }
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void ParametricClampedBar<T>
    ::AssembleNewMarkMatrix1()
    {
        Fill(T(0), Newmark_matrix_1_);
        for (int i = 0; i < Nx_; i++)
        {
            // Newmark's matrix at time n+1.
            Newmark_matrix_1_.Val(i, i) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(0, 0) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i, i) / Delta_t_ +
                0.5 * stiffness_FEM_matrix_(0, 0)
                * theta_stiffness_(theta_stiffness_index_(i));
            Newmark_matrix_1_.Val(i + 1, i + 1) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(1, 1) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i + 1, i + 1) / Delta_t_ +
                0.5 * stiffness_FEM_matrix_(1, 1)
                * theta_stiffness_(theta_stiffness_index_(i));
            Newmark_matrix_1_.Val(i, i + 1) +=
                theta_mass_(theta_mass_index_(i)) *
                2. * mass_FEM_matrix_(0, 1) / (Delta_t_ * Delta_t_) +
                damp_matrix_.Val(i, i + 1) / Delta_t_ +
                0.5 * stiffness_FEM_matrix_(0, 1)
                * theta_stiffness_(theta_stiffness_index_(i));
        }

        // Boundary condition by pseudo-elimination.
        Newmark_matrix_1_.Val(0, 0) = 1;
        Newmark_matrix_1_.Val(0, 1) = 0;
        Newmark_matrix_1_.Val(1, 0) = 0;

    }


    //! Assembles damp matrix.
    template <class T>
    void ParametricClampedBar<T>
    ::AssembleDampMatrix()
    {
        Fill(T(0), damp_matrix_);
        for (int i = 0; i < Nx_; i++)
        {
            damp_matrix_.Val(i, i) +=
                theta_damp_(theta_damp_index_(i)) *
                alpha_ * mass_FEM_matrix_(0, 0) +
                beta_ * stiffness_FEM_matrix_(0, 0)
                * theta_damp_(theta_damp_index_(i));
            damp_matrix_.Val(i + 1, i + 1) +=
                theta_damp_(theta_damp_index_(i)) *
                alpha_ * mass_FEM_matrix_(1, 1) +
                beta_ * stiffness_FEM_matrix_(1, 1)
                * theta_damp_(theta_damp_index_(i));
            damp_matrix_.Val(i, i + 1) +=
                theta_damp_(theta_damp_index_(i)) *
                alpha_ * mass_FEM_matrix_(0, 1) +
                beta_ * stiffness_FEM_matrix_(0, 1)
                * theta_damp_(theta_damp_index_(i));
        }

        // Boundary condition by pseudo-elimination.
        damp_matrix_.Val(0, 0) = 1;
        damp_matrix_.Val(0, 1) = 0;
        damp_matrix_.Val(1, 0) = 0;
    }


    //! Builds skeleton newmark, mass and damp matrices.
    template <class T>
    void ParametricClampedBar<T>
    ::AllocateSparseMatrix()
    {
        // Skeleton.
        int NvalSkel = 3 * (Nx_- 1) + 4;
        Vector<T> values_0(NvalSkel);
        values_0.Fill(T(0.));
        Vector<int> columns_0(NvalSkel);
        columns_0.Fill();
        Vector<int> rowindex_0(Ndof_ + 1);
        rowindex_0.Fill();
        columns_0(0) = 0;
        columns_0(1) = 1;
        rowindex_0(0) = 0;
        rowindex_0(1) = 2;

        for (int i = 1; i < Nx_; i++)
        {
            columns_0(3 * (i - 1) + 2) = i - 1;
            columns_0(3 * (i -1) + 3) = i;
            columns_0(3 * (i - 1) + 4) = i + 1;
            rowindex_0(i + 1) = rowindex_0(i) + 3;
        }

        columns_0(3 * (Nx_ - 1) + 2) = Ndof_ - 2;
        columns_0(3 * (Nx_ - 1) + 3) = Ndof_ - 1;
        rowindex_0(Ndof_) = 3 * (Nx_ - 1) + 4;

        // Store the upper part of the Newmark
        // matrix in a symmetric sparse data structure.
        int nnz = (NvalSkel + Ndof_) / 2;

        Vector<int> sym_col_0(nnz), sym_row_0(Ndof_ + 1);
        Vector<T> sym_values_0(nnz);

        int val_ind = 0, col_ind = 0;
        bool first_nz;
        for(int i = 0; i < Ndof_; i++)
        {
            first_nz = true;
            for(int j = rowindex_0(i); j < rowindex_0(i + 1); j++)
                if(columns_0(j) >= i)
                {
                    sym_values_0(val_ind) = values_0(j);
                    val_ind++;
                    sym_col_0(col_ind) = columns_0(j);
                    col_ind++;
                    if(first_nz)
                    {
                        sym_row_0(i)= col_ind - 1;
                        first_nz=false;
                    }
                }
        }

        sym_row_0(Ndof_) = nnz;

        Vector<T> sym_values_1(sym_values_0);
        Vector<int> sym_col_1(sym_col_0);
        Vector<int> sym_row_1(sym_row_0);
        Vector<T> sym_values_m(sym_values_0);
        Vector<int> sym_col_m(sym_col_0);
        Vector<int> sym_row_m(sym_row_0);
        Vector<T> sym_values_c(sym_values_0);
        Vector<int> sym_col_c(sym_col_0);
        Vector<int> sym_row_c(sym_row_0);

        // Remark: values, rowindex and colums are unlinked after used in
        // SetData therefore, we need one of each for each global matrix.
        Newmark_matrix_0_.SetData(Ndof_, Ndof_, sym_values_0,
                                  sym_row_0, sym_col_0);
        Newmark_matrix_1_.SetData(Ndof_, Ndof_, sym_values_1,
                                  sym_row_1, sym_col_1);
        mass_matrix_.SetData(Ndof_, Ndof_, sym_values_m,
                             sym_row_m, sym_col_m);
        damp_matrix_.SetData(Ndof_, Ndof_, sym_values_c,
                             sym_row_c, sym_col_c);
    }


}

#define VERDANDI_FILE_MODEL_PARAMETRICCLAMPEDBAR_CXX
#endif
