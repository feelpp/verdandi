// Copyright (C) 2008-2009 INRIA
// Author(s): Claire Mouton, Vivien Mallet
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_DIAGONALSPARSEOBSERVATIONMANAGER_CXX


#include <cstdlib>
#include "DiagonalSparseObservationManager.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! It entirely defines the operator: no dimension or size is associated
      with this implementation.
    */
    template <class T>
    DiagonalSparseObservationManager<T>
    ::DiagonalSparseObservationManager()
    {
    }


    //! Main constructor.
    /*! It defines the operator for a 2D regular grid.
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    DiagonalSparseObservationManager<T>
    ::DiagonalSparseObservationManager(const Model& model,
                                       string configuration_file):
        operator_sparse_(true), availability_(false), error_sparse_(true),
        error_matrix_availability_(true)
    {
        //   Initialize(model, configuration_file);
    }


    //! Destructor.
    template <class T>
    DiagonalSparseObservationManager<T>
    ::~DiagonalSparseObservationManager()
    {
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the observation manager.
    /*!
      param[in] model model.
      param[in] configuration_file configuration_file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    void DiagonalSparseObservationManager<T>
    ::Initialize(const Model& model, string configuration_file)
    {
        GetPot configuration_stream(configuration_file.c_str());

        //! First abscissa.
        double x_min_model = model.GetXMin();
        //! First ordinate.
        double y_min_model = model.GetYMin();

        //! Space step along x.
        double Delta_x_model = model.GetDeltaX();
        //! Space step along y.
        double Delta_y_model = model.GetDeltaY();

        Nx_model_ = model.GetNx();
        Ny_model_ = model.GetNy();

        configuration_stream.set_prefix("observation/");
        observation_file_ = configuration_stream("File",
                                                 "configuration_error");
        period_observation_ = configuration_stream("Period_observation", -1);
        Nskip_ = configuration_stream("Nskip", -1);
        error_variance_ = configuration_stream("Error_variance", -1.);

        Nobservation_ = Nx_model_ * Ny_model_;
        observation_.Reallocate(Nobservation_);

        /*** Building the sparse matrices ***/

        tangent_operator_matrix_.Initialize(Nobservation_, T(1));

        observation_error_covariance_matrix_.Initialize(Nobservation_,
                                                        error_variance_);
    }


    /////////////////////////////
    // OBSERVATIONS MANAGEMENT //
    /////////////////////////////


    //! Loads observations if they are available.
    /*!
      \param[in] model model.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    void DiagonalSparseObservationManager<T>
    ::LoadObservation(const Model& model)
    {
        int step = model.GetCurrentDate();

        availability_ = step % (period_observation_ * Nskip_) == 0;

        // If assimilation is finished (in prediction), no observation are
        // loaded.
        if (step > model.GetNtAssimilation())
            availability_ = false;

        if (availability_)
        {
            Matrix<double> input_data(Nx_model_, Ny_model_);
            ifstream file_stream;
            file_stream.open(observation_file_.c_str());

#ifdef SELDON_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw IOError("DiagonalSparseObservationManager"
                              "::LoadObservation(model)",
                              string("Unable to open file \"")
                              + observation_file_ + "\".");
#endif

            streampos position = step / period_observation_
                * Nx_model_ * Ny_model_ * sizeof(double);
            file_stream.seekg(position);

            // To be optimized: use a method reading only one step instead.
            input_data.Read(file_stream, false);

            file_stream.close();

            // To be optimized: avoid the copy by reading directly to a data
            // container.
            for (int i = 0; i < Nobservation_; i++)
            {
                div_t division;
                division = div(i, Ny_model_);
                observation_(i) = input_data(division.quot, division.rem);
            }
        }
    }


    //! Checks whether observations are available.
    /*!
      \return True if observations are available, false otherwise.
    */
    template <class T>
    bool DiagonalSparseObservationManager<T>::HasObservation() const
    {
        return availability_;
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the number of observations.
    /*!
      \return Observation dimension.
    */
    template <class T>
    int DiagonalSparseObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
    }


    //! Returns the values of observations.
    /*!
      \return The values of observations.
    */
    template <class T>
    const Vector<T>& DiagonalSparseObservationManager<T>
    ::GetObservation() const
    {
        return observation_;
    }


    /*! \brief Checks whether the observation operator is available in a
      sparse matrix. */
    /*!
      \return True if the observation operator is available in a sparse
      matrix, false otherwise.
    */
    template <class T>
    bool DiagonalSparseObservationManager<T>::IsOperatorSparse() const
    {
        return operator_sparse_;
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool DiagonalSparseObservationManager<T>::IsErrorSparse() const
    {
        return error_sparse_;
    }


    /*! \brief Checks whether the observation error covariance is available in
      a matrix. */
    /*!
      \return True if the observation error covariance is available in a
      matrix, false otherwise.
    */
    template <class T>
    bool DiagonalSparseObservationManager<T>::HasErrorMatrix() const
    {
        return error_matrix_availability_;
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::ApplyOperator(const Vector<T>& x, Vector<T>& y) const
    {
        y = x;
    }


    //! Applies the tangent linear operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::ApplyTangentOperator(const Vector<T>& x, Vector<T>& y) const
    {
        ApplyOperator(x, y);
    }


    //! Linearized observation operator.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the linearized operator.
    */
    template <class T>
    T DiagonalSparseObservationManager<T>
    ::GetTangentOperator(int i, int j) const
    {
        return tangent_operator_matrix_.GetValue(i, j);
    }


    //! Linearized observation operator.
    /*!
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the linearized
      operator.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::GetTangentOperatorRow(int row, Vector<T>& tangent_operator_row) const
    {
        tangent_operator_matrix_.GetRow(row, tangent_operator_row);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& DiagonalSparseObservationManager<T>
    ::GetTangentOperatorMatrix() const
    {
        return tangent_operator_matrix_.GetMatrix();
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::ApplyAdjointOperator(const Vector<T>& x, Vector<T>& y) const
    {
        y = x;
    }


    //! Gets innovation.
    /*!
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::GetInnovation(const Vector<T>& state, Vector<T>& innovation) const
    {
        ApplyOperator(state, innovation);
        Mlt(T(-1), innovation);
        Add(T(1), GetObservation(), innovation);
    }


    //! Checks whether a BLUE correction is available.
    /*!
      \return True if a BLUE correction is available, false otherwise.
    */
    template <class T>
    bool DiagonalSparseObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined("DiagonalSparseObservationManager"
                             "::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void DiagonalSparseObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined("DiagonalSparseObservationManager"
                             "::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T DiagonalSparseObservationManager<T>
    ::GetObservationErrorCovariance(int i, int j) const
    {
        return observation_error_covariance_matrix_.GetValue(i, j);
    }


    //! Observation error covariance matrix.
    /*!
      \return The matrix of the observation error covariance.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& DiagonalSparseObservationManager<T>
    ::GetObservationErrorCovarianceMatrix() const
    {
        return observation_error_covariance_matrix_.GetMatrix();
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_DIAGONALSPARSEOBSERVATIONMANAGER_CXX
#endif
