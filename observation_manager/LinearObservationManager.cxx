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
        availability_(false)
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
        configuration_stream.set("File", observation_file_);
        configuration_stream.set("Period_observation",
                                 period_observation_,  "> 0");
        configuration_stream.set("Nskip", Nskip_, "> 0");

        configuration_stream.set("error/Variance", error_variance_value_,
                                 "> 0");

        configuration_stream.set("operator/Definition",
                                 operator_definition_, "'diagonal'| 'file'");
        configuration_stream.set("operator/Diagonal_value",
                                 operator_diagonal_value_);
        configuration_stream.set("operator/File", operator_file_);

        Nobservation_ = Nx_model_ * Ny_model_;
        observation_.Reallocate(Nobservation_);

        /*** Building the matrices ***/

#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
        build_diagonal_sparse_matrix(Nobservation_, operator_diagonal_value_,
                                     tangent_operator_matrix_);
#endif

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
                                     error_variance_);
#endif

#ifdef VERDANDI_TANGENT_OPERATOR_DENSE
        if (strcmp(operator_definition_.c_str(), "file") == 0)
            tangent_operator_matrix_.Read(operator_file_);
#endif
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
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
        return true;
#else
        return false;
#endif
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>::IsErrorSparse() const
    {
#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        return true;
#else
        return false;
#endif
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
#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        return true;
#else
        return false;
#endif
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
    template <class state_vector>
    void LinearObservationManager<T>
    ::ApplyOperator(const state_vector& x, Vector<T>& y) const
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
    template <class state_vector>
    void LinearObservationManager<T>
    ::ApplyTangentOperator(const state_vector& x, Vector<T>& y) const
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
    ::GetTangentOperatorRow(int row, typename LinearObservationManager<T>
                            ::tangent_operator_row& tangent_operator_row)
        const
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
    const typename LinearObservationManager<T>
    ::tangent_operator_matrix& LinearObservationManager<T>
    ::GetTangentOperatorMatrix() const
    {
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
        if (strcmp(operator_definition_.c_str(), "diagonal") == 0)
            return tangent_operator_matrix_;

        // Dense operator or operator defined in a file.
        else
#endif
            throw ErrorUndefined("LinearObservationManager"
                                 "::GetTangentOperatorMatrix()");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    template <class state_vector>
    void LinearObservationManager<T>
    ::ApplyAdjointOperator(const state_vector& x, Vector<T>& y) const
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
    template <class state_vector>
    void LinearObservationManager<T>
    ::GetInnovation(const state_vector& state, Vector<T>& innovation) const
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
            return error_variance_value_;
        else
            return T(0);
    }


    //! Observation error covariance matrix.
    /*!
      \return The matrix of the observation error covariance.
    */
    template <class T>
    const typename LinearObservationManager<T>
    ::error_variance& LinearObservationManager<T>
    ::GetObservationErrorVariance() const
    {
#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        return error_variance_;
#else
        throw ErrorUndefined("LinearObservationManager"
                             "::GetObservationErrorVariance()");
#endif
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string LinearObservationManager<T>::GetName() const
    {
        return "LinearObservationManager";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void LinearObservationManager<T>::Message(string message)
    {
        cout << "OM received: " << message << endl;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_LINEAROBSERVATIONMANAGER_CXX
#endif
