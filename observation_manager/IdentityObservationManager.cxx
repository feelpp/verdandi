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


#ifndef VERDANDI_FILE_IDENTITYOBSERVATIONMANAGER_CXX

#include <cstdlib>
#include "IdentityObservationManager.hxx"


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
    IdentityObservationManager<T>::IdentityObservationManager()
    {
    }


    //! Main constructor.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam ClassModel the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class ClassModel>
    IdentityObservationManager<T>
    ::IdentityObservationManager(const ClassModel& model,
                                 string configuration_file):
        operator_sparse_(false), availability_(false), error_sparse_(false),
        error_matrix_availability_(false)
    {
        //   Initialize(model, configuration_file);
    }


    //! Destructor.
    template <class T>
    IdentityObservationManager<T>
    ::~IdentityObservationManager()
    {
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the observation manager.
    /*!
      param[in] model model.
      param[in] configuration_file configuration_file.
      \tparam ClassModel the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class ClassModel>
    void IdentityObservationManager<T>
    ::Initialize(const ClassModel& model, string configuration_file)
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
    }


    /////////////////////////////
    // OBSERVATIONS MANAGEMENT //
    /////////////////////////////


    //! Loads observations if they are available.
    /*!
      \param[in] model model.
      \tparam ClassModel the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class ClassModel>
    void IdentityObservationManager<T>
    ::LoadObservation(const ClassModel& model)
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
                throw IOError("IdentityObservationManager"
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
    bool IdentityObservationManager<T>::HasObservation() const
    {
        return availability_;
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the observation number.
    /*!
      \return Observation dimension.
    */
    template <class T>
    int IdentityObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
    }


    //! Returns the values of observations.
    /*!
      \return The values of observations.
    */
    template <class T>
    const Vector<T>& IdentityObservationManager<T>
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
    bool IdentityObservationManager<T>::IsOperatorSparse() const
    {
        return operator_sparse_;
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool IdentityObservationManager<T>::IsErrorSparse() const
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
    bool IdentityObservationManager<T>::HasErrorMatrix() const
    {
        return error_matrix_availability_;
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the operator (thus, the identity) to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    void IdentityObservationManager<T>
    ::ApplyOperator(const Vector<T>& x, Vector<T>& y) const
    {
        y = x;
    }


    /*! \brief Applies the tangent linear operator (thus, the identity) to a
      given vector. */
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x.
    */
    template <class T>
    void IdentityObservationManager<T>
    ::ApplyTangentOperator(const Vector<T>& x, Vector<T>& y) const
    {
        ApplyOperator(x, y);
    }


    //! Linearized observation operator, which is simply the identity.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the linearized operator: 1 if \a i
      equals \a j, and 0 otherwise.
    */
    template <class T>
    T IdentityObservationManager<T>
    ::GetTangentOperator(int i, int j) const
    {
        if (i == j)
            return T(1);
        else
            return T(0);
    }


    //! Linearized observation operator.
    /*!
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the linearized
      operator.
    */
    template <class T>
    void IdentityObservationManager<T>
    ::GetTangentOperatorRow(int row, Vector<T>& tangent_operator_row) const
    {
        tangent_operator_row.Reallocate(Nobservation_);
        tangent_operator_row.Zero();
        tangent_operator_row(row) = T(1);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& IdentityObservationManager<T>
    ::GetTangentOperatorMatrix() const
    {
        throw ErrorUndefined("IdentityObservationManager"
                             "::GetTangentOperatorMatrix()");
    }


    //! Applies the adjoint operator (thus, the identity) to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    void IdentityObservationManager<T>
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
    void IdentityObservationManager<T>
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
    bool IdentityObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined("IdentityObservationManager"
                             "::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void IdentityObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined("IdentityObservationManager"
                             "::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T IdentityObservationManager<T>
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
    const Matrix<T>& IdentityObservationManager<T>
    ::GetObservationErrorCovarianceMatrix() const
    {
        throw ErrorUndefined("IdentityObservationManager"
                             "::GetObservationErrorCovarianceMatrix()");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_CXX
#endif
