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


#ifndef VERDANDI_FILE_LINEAROBSERVATIONMANAGER_CXX

#include <cstdlib>
#include "LinearObservationManager.hxx"


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
    LinearObservationManager<T>::LinearObservationManager()
    {
    }


    //! Main constructor.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    LinearObservationManager<T>
    ::LinearObservationManager(const Model& model,
                               string configuration_file):
        operator_dense_(false), availability_(false)
    {
        //   Initialize(model, configuration_file);
    }


    //! Destructor.
    template <class T>
    LinearObservationManager<T>
    ::~LinearObservationManager()
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
    void LinearObservationManager<T>
    ::Initialize(const Model& model, string configuration_file)
    {
        GetPot configuration_stream(configuration_file);

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
        configuration_stream.put("File", observation_file_);
        configuration_stream.put("Period_observation",  "> 0",
                                 period_observation_);
        configuration_stream.put("Nskip", "> 0", Nskip_);

        configuration_stream.put("error/Variance", "> 0", error_variance_);
        configuration_stream.put("error/Sparse", error_sparse_);
        configuration_stream.put("error/Matrix_availability",
                                 error_matrix_availability_);

        configuration_stream.put("operator/Sparse", operator_sparse_);
        configuration_stream.put("operator/Definition", "'diagonal'| file'",
                                 operator_definition_);
        configuration_stream.put("operator/Diagonal_value",
                                 operator_diagonal_value_);
        configuration_stream.put("operator/File", operator_file_);

        Nobservation_ = Nx_model_ * Ny_model_;
        observation_.Reallocate(Nobservation_);

        /*** Building the matrices ***/

        if (operator_sparse_)
            tangent_operator_sparse_matrix_
                .Initialize(Nobservation_, operator_diagonal_value_);

        if (error_sparse_)
            observation_error_covariance_matrix_.Initialize(Nobservation_,
                                                            error_variance_);

        if (strcmp(operator_definition_.c_str(), "file") == 0)
        {
            operator_dense_ = true;
            tangent_operator_matrix_.Read(operator_file_);
        }
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
    void LinearObservationManager<T>
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
                throw IOError("LinearObservationManager"
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
    bool LinearObservationManager<T>::HasObservation() const
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
    int LinearObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
    }


    //! Returns the values of observations.
    /*!
      \return The values of observations.
    */
    template <class T>
    const Vector<T>& LinearObservationManager<T>
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
    bool LinearObservationManager<T>::IsOperatorSparse() const
    {
        return operator_sparse_;
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>::IsErrorSparse() const
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
    bool LinearObservationManager<T>::HasErrorMatrix() const
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
    void LinearObservationManager<T>
    ::ApplyOperator(const Vector<T>& x, Vector<T>& y) const
    {
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0)
        {
            y = x;
            Mlt(operator_diagonal_value_, y);
        }

        // Operator defined in a file.
        else
            Mlt(tangent_operator_matrix_, x, y);
    }


    //! Applies the tangent linear operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x.
    */
    template <class T>
    void LinearObservationManager<T>
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
    T LinearObservationManager<T>
    ::GetTangentOperator(int i, int j) const
    {
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0)
        {
            if (i == j)
                return operator_diagonal_value_;
            else
                return T(0);
        }

        // Operator defined in a file.
        else
            return tangent_operator_matrix_(i, j);
    }


    //! Linearized observation operator.
    /*!
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the linearized
      operator.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetTangentOperatorRow(int row, Vector<T>& tangent_operator_row) const
    {
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0)
        {
            // if (operator_sparse_)
            // {
            //     tangent_operator_row.Reallocate(1);
            //     tangent_operator_row.Index(0) = row;
            //     tangent_operator_row.Fill(operator_diagonal_value_);
            // }

            // // Dense operator.
            // else
            {
                tangent_operator_row.Reallocate(Nobservation_);
                tangent_operator_row.Zero();
                tangent_operator_row(row) = operator_diagonal_value_;
            }
        }

        // Operator defined in a file.
        else
            GetRow(tangent_operator_matrix_, row, tangent_operator_row);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& LinearObservationManager<T>
    ::GetTangentOperatorMatrix() const
    {
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0 and
            operator_sparse_)
            return tangent_operator_sparse_matrix_.GetMatrix();

        // Dense operator or operator defined in a file.
        else
            throw ErrorUndefined("LinearObservationManager"
                                 "::GetTangentOperatorMatrix()");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ApplyAdjointOperator(const Vector<T>& x, Vector<T>& y) const
    {
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0)
        {
            y = x;
            Mlt(operator_diagonal_value_, y);
        }

        // Operator defined in a file.
        else
            Mlt(1., SeldonTrans, tangent_operator_matrix_, x, 0., y);
    }


    //! Gets innovation.
    /*!
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class T>
    void LinearObservationManager<T>
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
    bool LinearObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined("LinearObservationManager"
                             "::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined("LinearObservationManager"
                             "::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T LinearObservationManager<T>
    ::GetObservationErrorCovariance(int i, int j) const
    {
        if (i == j)
            return error_variance_;
        else
            return T(0);
    }


    //! Observation error covariance matrix.
    /*!
      \return The matrix of the observation error covariance.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& LinearObservationManager<T>
    ::GetObservationErrorCovarianceMatrix() const
    {
        if (error_matrix_availability_)
            return observation_error_covariance_matrix_.GetMatrix();
        else
            throw ErrorUndefined("LinearObservationManager"
                                 "::GetObservationErrorCovarianceMatrix()");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_LINEAROBSERVATIONMANAGER_CXX
#endif
