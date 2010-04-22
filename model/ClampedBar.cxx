// Copyright (C) 2008-2009 INRIA
// Author(s): Dominique Chapelle, Philippe Moireau
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


#ifndef VERDANDI_FILE_MODEL_CLAMPEDBAR_CXX


#include "ClampedBar.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////

    //! Constructor.
    template <class T>
    ClampedBar<T>::ClampedBar()
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    ClampedBar<T>::ClampedBar(string configuration_file)
    {
    }


    //! Destructor.
    template <class T>
    ClampedBar<T>::~ClampedBar()
    {
    }

    ////////////////
    // INITIALIZE //
    ////////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ClampedBar<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        GetPot configuration_stream(configuration_file, "#", "\n");

        configuration_stream.set_prefix("clamped_bar/domain/");
        configuration_stream.set("Bar_length", bar_length_);
        configuration_stream.set("Nx", Nx_);
        configuration_stream.set("Delta_t", Delta_t_);
        configuration_stream.set("Final_date", final_date_);

        configuration_stream.set_prefix("clamped_bar/error_statistics/");
        configuration_stream.set("Background_error_variance",
                                 background_error_variance_value_, ">= 0");
        configuration_stream.set("Background_error_scale",
                                 Balgovind_scale_background_, "> 0");
        
        Ndof_ = Nx_ + 1;

#ifdef VERDANDI_BACKGROUND_ERROR_SPARSE
        build_diagonal_sparse_matrix(GetNstate(),
                                     background_error_variance_value_,
                                     background_error_variance_);
#endif

        configuration_stream.set_prefix("clamped_bar/physics/");
        configuration_stream.set("Young_modulus", Young_modulus_);
        configuration_stream.set("mass_density", mass_density_);

        /*** Ouput saver ***/

        output_saver_.Initialize(configuration_file,
                                 "clamped_bar/output_saver/");
        output_saver_.Empty("disp_0");
        output_saver_.Empty("velo_0");

        /*** Allocation ***/

        Delta_x_ = bar_length_ / Nx_;
        date_vector_.reserve(floor(final_date_ / Delta_t_));
        date_vector_ = vector<double>(1, 0.);
        date_ = 0.;

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

        Logger::StdOut("values_0", values_0);
        Logger::StdOut("columns_0", columns_0);
        Logger::StdOut("rowindex_0", rowindex_0);

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

        //Remark: values,rowindex and colums are unlinked after used in
        // SetData therefore, we need one of each for each global matrice.
        Newmark_matrix_0_.SetData(Ndof_, Ndof_, sym_values_0,
                                  sym_row_0, sym_col_0);
        Newmark_matrix_1_.SetData(Ndof_, Ndof_, sym_values_1,
                                  sym_row_1, sym_col_1);
        Mass_matrix_.SetData(Ndof_, Ndof_, sym_values_m,
                             sym_row_m, sym_col_m);

        // Creates local matrix.
        Mass_matrix_el_.Reallocate(2, 2);
        Stiff_matrix_el_.Reallocate(2, 2);

        InitializeFirstStep();
    }


    //! Initializes the first time step for the model.
    template <class T>
    void ClampedBar<T>::InitializeFirstStep()
    {
        /*** Build model ***/
        // Elementary mass matrix construction.
        T mass_lin = T(mass_density_ * Delta_x_ / 3);
        Mass_matrix_el_(0, 0) = mass_lin;
        Mass_matrix_el_(1, 1) = mass_lin;
        Mass_matrix_el_(0, 1) = mass_lin / 2;

        // Elementary stifness matrix construction.
        T stiff_lin = T(Young_modulus_ / Delta_x_);
        Stiff_matrix_el_(0, 0) = stiff_lin;
        Stiff_matrix_el_(1, 1) = stiff_lin;
        Stiff_matrix_el_(0, 1) = -stiff_lin;

        // Assembling process.
        for (int i=0; i < Nx_; i++)
        {
            //Mass Matrice.
            Mass_matrix_.Val(i, i) += Mass_matrix_el_(0, 0);
            Mass_matrix_.Val(i + 1, i + 1) += Mass_matrix_el_(1, 1);
            Mass_matrix_.Val(i, i + 1) += Mass_matrix_el_(0, 1);

            // Newmark's matrice at time n.
            Newmark_matrix_0_.Val(i, i) +=
                2. * Mass_matrix_el_(0, 0) / (Delta_t_ * Delta_t_) +
                -0.5 * Stiff_matrix_el_(0, 0);
            Newmark_matrix_0_.Val(i + 1, i + 1) +=
                2. * Mass_matrix_el_(1, 1) / (Delta_t_ * Delta_t_) +
                -0.5 * Stiff_matrix_el_(1, 1);
            Newmark_matrix_0_.Val(i, i + 1) +=
                2. * Mass_matrix_el_(0, 1) / (Delta_t_ * Delta_t_) +
                -0.5 * Stiff_matrix_el_(0, 1);

            // Newmark's matrice at time n+1.
            Newmark_matrix_1_.Val(i, i) +=
                2. * Mass_matrix_el_(0, 0) / (Delta_t_ * Delta_t_) +
                0.5 * Stiff_matrix_el_(0, 0);
            Newmark_matrix_1_.Val(i + 1, i + 1) +=
                2. * Mass_matrix_el_(1, 1) / (Delta_t_ * Delta_t_) +
                0.5 * Stiff_matrix_el_(1, 1);
            Newmark_matrix_1_.Val(i, i + 1) +=
                2. * Mass_matrix_el_(0, 1) / (Delta_t_ * Delta_t_) +
                0.5 * Stiff_matrix_el_(0, 1);

        }


        Logger::StdOut(*this, "Assembled matrices");
        Logger::StdOut("Mass_matrix_", Mass_matrix_);
        Logger::StdOut("Newmark_matrix_1_", Newmark_matrix_1_);
        Logger::StdOut("Newmark_matrix_0_", Newmark_matrix_0_);

        // Dirichlet conditions (should be better to read in the skeleton).
        Newmark_matrix_1_.Val(0, 0) = 1;
        Newmark_matrix_1_.Val(0, 1) = 0;
        Newmark_matrix_1_.Val(1, 0) = 0;

        // Vector initialization.
        disp_0_.Reallocate(Ndof_);
        disp_1_.Reallocate(Ndof_);
        velo_0_.Reallocate(Ndof_);
        velo_1_.Reallocate(Ndof_);
        force_.Reallocate(Ndof_);
        disp_0_.Fill(T(0.));
        velo_0_.Fill(T(0.));
        disp_1_.Fill(T(0.));
        velo_1_.Fill(T(0.));
        force_.Fill(T(0.));

        // Initial condition.
        for (int i = 0; i < Ndof_; i++)
            disp_0_(i) = T(i) / T(Ndof_ - 1);

#ifdef VERDANDI_WITH_DIRECT_SOLVER
        GetLU(Newmark_matrix_1_, mat_lu, true);
#endif

        Logger::StdOut("disp_0_", disp_0_);
        Logger::StdOut("velo_0_", velo_0_);
        Logger::StdOut("force_", force_);
    }


    //! Initializes the current time step for the model.
    template <class T>
    void ClampedBar<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    template <class T>
    void ClampedBar<T>::Forward()
    {
        // Update time.
        date_ += Delta_t_;
        date_vector_.push_back(date_);

        Logger::StdOutCommand("hline", "=");
        Logger::StdOut("Time", date_);
        Logger::StdOutCommand("hline", "-");

        // Right hand side.
        force_.Fill(T(0.));

        MltAdd(2. / Delta_t_, Mass_matrix_, velo_0_, 1, force_);
        MltAdd(1., Newmark_matrix_0_, disp_0_, 1, force_);

        // Dirichlet conditions .
        force_(0) = 0;

#ifdef VERDANDI_WITH_DIRECT_SOLVER
        Vector<T> tmp(Ndof_);
        tmp = force_;
        SolveLU(mat_lu, tmp);
        disp_1_ = tmp;
#else
        // Initialization of the Gmres parameters.
        int nb_max_iter = 1000;
        double tolerance = 1e-6;
        Iteration<double> iter(nb_max_iter, tolerance);
        Preconditioner_Base precond;
        // No preconditioning.
        iter.SetRestart(5);
        Gmres(Newmark_matrix_1_, disp_1_, force_, precond, iter);
#endif

        velo_1_.Fill(T(0.));
        Add(-1., velo_0_, velo_1_);
        Add(2. / Delta_t_, disp_1_, velo_1_);
        Add(-2. / Delta_t_, disp_0_, velo_1_);

        // Update.
        disp_0_ = disp_1_;
        velo_0_ = velo_1_;

        Logger::StdOut("disp_0_", disp_0_);
        Logger::StdOut("velo_0_", velo_0_);
        Logger::StdOut("force_", force_);
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool ClampedBar<T>::HasFinished() const
    {
        return date_ >= final_date_;
    }


    //! Saves the simulated data.
    /*! It saves the displacement 'disp_0_' and  the velocity 'velo_0_'.
     */
    template <class T>
    void ClampedBar<T>::Save()
    {
        output_saver_.Save(disp_0_, date_, "disp_0");
        output_saver_.Save(velo_0_, date_, "velo_0");
    }


    //! Returns the current date.
    /*!
      \return The current date.
    */
    template <class T>
    double ClampedBar<T>::GetDate() const
    {
        return date_;
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    template <class T>
    int ClampedBar<T>::GetNstate() const
    {
        return 2 * Ndof_;
    }


    //! Provides the reduced state vector.
    /*!
      \param[out] state the reduced state vector.
    */
    template <class T>
    void ClampedBar<T>::GetState(state_vector& state) const
    {
        int position = 0;
        state.Reallocate(2 * Ndof_);
        for (int i = 0; i < Ndof_; i++)
            state(position++) = disp_0_(i);
        for (int i = 0; i < Ndof_; i++)
            state(position++) = velo_0_(i);
    }


    //! Sets the reduced state vector.
    /*! Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
    */
    template <class T>
    void ClampedBar<T>::SetState(state_vector& state)
    {
        int position = 0;
        for (int i = 0; i < Ndof_; i++)
            disp_0_(i) = state(position++);
        for (int i = 0; i < Ndof_; i++)
            velo_0_(i) = state(position++);
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    void ClampedBar<T>::GetFullState(state_vector& state) const
    {
        throw ErrorUndefined("ClampedBar"
                             "::GetFullState()");
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    template <class T>
    void ClampedBar<T>::SetFullState(const state_vector& state)
    {
        throw ErrorUndefined("ClampedBar"
                             "::SetFullState()");
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] error_covariance_row the value of row number \a row.
    */
    template <class T>
    void ClampedBar<T>
    ::GetBackgroundErrorCovarianceRow(int row, error_covariance_row&
                                      error_covariance_row)
    {
#ifdef VERDANDI_BACKGROUND_ERROR_SPARSE
        {
            error_covariance_row.Reallocate(GetNstate());
            error_covariance_row.Zero();
            error_covariance_row(row) = background_error_variance_value_;
        }
# else
        {
#ifdef VERDANDI_BACKGROUND_ERROR_DENSE
            {
                error_covariance_row.Reallocate(GetNstate());
                error_covariance_row.Zero();
                error_covariance_row(row)
                    = background_error_variance_value_;
            }
#else
            {
                // The row has already been computed.
                if (row == current_row_)
                    error_covariance_row = error_covariance_row_;
                else
                {
                    int i;
                    current_row_ = row;
                    current_column_ = -1;

                    // Positions related to 'row'.
                    int i_row = row;
                    int j_row = row - i_row;

                    T distance_x;
                    T distance;
                    int position = 0;
                    for (i = 0; i < Nx_; i++)
                    {
                        distance_x = Delta_x_ * T(i - i_row);
                        distance = sqrt(distance_x * distance_x)
                            / Balgovind_scale_background_;
                        error_covariance_row_(position++)
                            = background_error_variance_value_
                            * (1. + distance) * exp(-distance);
                    }
                    error_covariance_row = error_covariance_row_;
                }
            }
#endif
        }
#endif
    }


    //! Returns the background error covariance matrix (B) if available.
    /*! Returns the background error covariance matrix (B) if available,
      raises an exception otherwise.
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename ClampedBar<T>::background_error_variance& ClampedBar<T>
    ::GetBackgroundErrorVarianceMatrix() const
    {
#ifdef VERDANDI_BACKGROUND_ERROR_SPARSE
        return background_error_variance_;
#else
        throw ErrorUndefined(
            "ClampedBar::GetBackgroundErrorVarianceMatrix()",
            "the background error covariance matrix is not available!");
#endif
    }


    //! Checks if the error covariance matrix is sparse.
    /*!
      \return True if there is a sparse error matrix, false otherwise.
    */
    template <class T>
    bool ClampedBar<T>::IsErrorSparse() const
    {
#ifdef VERDANDI_BACKGROUND_ERROR_SPARSE
        return true;
#else
        return false;
#endif

    }


    //! Returns the name of the class.
    template <class T>
    string ClampedBar<T>::GetName() const
    {
        return "ClampedBar";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void ClampedBar<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


}

#define VERDANDI_FILE_MODEL_CLAMPEDBAR_CXX
#endif
