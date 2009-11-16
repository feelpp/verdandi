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


#ifndef VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_CXX


#include "GridToNetworkObservationManager.hxx"


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
    GridToNetworkObservationManager<T>
    ::GridToNetworkObservationManager()
    {
    }


    //! Main constructor.
    /*! It defines the operator for a 2D regular grid.
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam ClassModel the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class ClassModel>
    GridToNetworkObservationManager<T>
    ::GridToNetworkObservationManager(const ClassModel& model,
                                      string configuration_file):
        availability_(false)
    {
        //   Initialize(model, configuration_file);
    }


    //! Destructor.
    template <class T>
    GridToNetworkObservationManager<T>
    ::~GridToNetworkObservationManager()
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
    void GridToNetworkObservationManager<T>
    ::Initialize(const ClassModel& model, string configuration_file)
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
                                 period_observation_, "> 0");
        configuration_stream.set("Nskip", Nskip_, "> 0");
        configuration_stream.set("error/Variance", error_variance_value_,
                                 "> 0");

        configuration_stream.set_prefix("observation/location/");

        string observation_location;
        configuration_stream.set("Observation_location",
                                 observation_location);
        vector<string> observation_location_vector
            = split(observation_location);
        int value;
        for (int i = 0; i < int(observation_location_vector.size() - 1);
             i += 2)
        {
            to_num(observation_location_vector[i], value);
            location_x_.PushBack(value);
            to_num(observation_location_vector[i + 1], value);
            location_y_.PushBack(value);
        }

        Nobservation_ = int(location_x_.GetSize());

        for (int i = 0; i < Nobservation_; i++)
            if (location_x_(i) < 0 || location_x_(i) >= Nx_model_)
                throw "Index along x should be in [0, "
                    + to_str(Nx_model_ - 1)
                    + "], but " + to_str(location_x_(i)) + " was provided.";
        for (int i = 0; i < Nobservation_; i++)
            if (location_y_(i) < 0 || location_y_(i) >= Ny_model_)
                throw "Index along y should be in [0, "
                    + to_str(Ny_model_ - 1)
                    + "], but " + to_str(location_y_(i)) + " was provided.";

        observation_.Reallocate(Nobservation_);

        interpolation_index_.Reallocate(Nobservation_, 4);
        interpolation_weight_.Reallocate(Nobservation_, 4);

        int i, j;
        T weight_x, weight_y;
        for (int location = 0; location < Nobservation_; location++)
        {
            i = int((location_x_(location) - x_min_model) / Delta_x_model);
            j = int((location_y_(location) - y_min_model) / Delta_y_model);
            interpolation_index_(location, 0) = i * Ny_model_ + j;
            interpolation_index_(location, 1) = i * Ny_model_ + j + 1;
            interpolation_index_(location, 2) = (i + 1) * Ny_model_ + j;
            interpolation_index_(location, 3) = (i + 1) * Ny_model_ + j + 1;

            weight_x = (location_x_(location) - x_min_model
                        - T(i) * Delta_x_model) / Delta_x_model;
            weight_y = (location_y_(location) - y_min_model
                        - T(j) * Delta_y_model) / Delta_y_model;
            interpolation_weight_(location, 0) = (T(1) - weight_x)
                * (T(1) - weight_y);
            interpolation_weight_(location, 1) = (T(1) - weight_x) * weight_y;
            interpolation_weight_(location, 2) = weight_x * (T(1) - weight_y);
            interpolation_weight_(location, 3) = weight_x * weight_y;
        }

        // Sets all active by default.
        SetAllActive();
    }


    //! Sets all observation locations as active.
    template <class T>
    void GridToNetworkObservationManager<T>::SetAllActive()
    {
        active_interpolation_index_ = interpolation_index_;
        active_interpolation_weight_ = interpolation_weight_;
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
    void GridToNetworkObservationManager<T>
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
                throw IOError("GridToNetworkObservationManager" \
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
                observation_(i) = input_data(location_x_(i), location_y_(i));
        }
    }


    //! Checks whether observations are available.
    /*!
      \return True if observations are available, false otherwise.
    */
    template <class T>
    bool GridToNetworkObservationManager<T>::HasObservation() const
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
    int GridToNetworkObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
    }


    //! Returns the values of observations.
    /*!
      \return The values of observations.
    */
    template <class T>
    const Vector<T>& GridToNetworkObservationManager<T>
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
    bool GridToNetworkObservationManager<T>::IsOperatorSparse() const
    {
        return false;
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool GridToNetworkObservationManager<T>::IsErrorSparse() const
    {
        return false;
    }


    /*! \brief Checks whether the observation error covariance is available in
      a matrix. */
    /*!
      \return True if the observation error covariance is available in a
      matrix, false otherwise.
    */
    template <class T>
    bool GridToNetworkObservationManager<T>::HasErrorMatrix() const
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
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class T>
    template <class state_vector>
    void GridToNetworkObservationManager<T>
    ::ApplyOperator(const state_vector& x, Vector<T>& y) const
    {
        y.Reallocate(active_interpolation_index_.GetM());
        y.Zero();
        for (int i = 0; i < y.GetLength(); i++)
            y(i) = active_interpolation_weight_(i, 0)
                * x(active_interpolation_index_(i, 0))
                + active_interpolation_weight_(i, 1)
                * x(active_interpolation_index_(i, 1))
                + active_interpolation_weight_(i, 2)
                * x(active_interpolation_index_(i, 2))
                + active_interpolation_weight_(i, 3)
                * x(active_interpolation_index_(i, 3));
    }


    //! Applies the tangent linear operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x. It is
      resized if needed.
    */
    template <class T>
    template <class state_vector>
    void GridToNetworkObservationManager<T>
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
    T GridToNetworkObservationManager<T>
    ::GetTangentOperator(int i, int j) const
    {
        if (j == active_interpolation_index_(i, 0))
            return active_interpolation_weight_(i, 0);
        else if (j == active_interpolation_index_(i, 1))
            return active_interpolation_weight_(i, 1);
        else if (j == active_interpolation_index_(i, 2))
            return active_interpolation_weight_(i, 2);
        else if (j == active_interpolation_index_(i, 3))
            return active_interpolation_weight_(i, 3);
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
    void GridToNetworkObservationManager<T>
    ::GetTangentOperatorRow(int row, GridToNetworkObservationManager<T>
                            ::tangent_operator_row& tangent_operator_row)
        const
    {
        for (int i = 0; i < tangent_operator_row.GetLength(); i++)
            tangent_operator_row(i) = GetTangentOperator(row, i);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const typename
    GridToNetworkObservationManager<T>::tangent_operator_matrix&
    GridToNetworkObservationManager<T>
    ::GetTangentOperatorMatrix() const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::GetTangentOperatorMatrix()");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class T>
    template <class state_vector>
    void GridToNetworkObservationManager<T>
    ::ApplyAdjointOperator(const state_vector& x, Vector<T>& y) const
    {
        y.Reallocate(Nx_model_ * Ny_model_);
        y.Zero();
        for (int i = 0; i < Nx_model_ * Ny_model_; i++)
            for (int j = 0; j < active_interpolation_index_.GetM(); j++)
                if (i == active_interpolation_index_(j, 0))
                    y(i) += active_interpolation_weight_(j, 0) * x(j);
                else if (i == active_interpolation_index_(j, 1))
                    y(i) += active_interpolation_weight_(j, 1) * x(j);
                else if (i == active_interpolation_index_(j, 2))
                    y(i) += active_interpolation_weight_(j, 2) * x(j);
                else if (i == active_interpolation_index_(j, 3))
                    y(i) += active_interpolation_weight_(j, 3) * x(j);
    }


    //! Gets innovation.
    /*!
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class T>
    template <class state_vector>
    void GridToNetworkObservationManager<T>
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
    bool GridToNetworkObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T GridToNetworkObservationManager<T>
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
    const typename
    GridToNetworkObservationManager<T>::error_variance&
    GridToNetworkObservationManager<T>
    ::GetObservationErrorVariance() const
    {
        throw ErrorUndefined(
            string("GridToNetworkObservationManager")
            +"::GetObservationErrorVariance()");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_CXX
#endif
