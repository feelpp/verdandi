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


#ifndef VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_HXX

#include <iostream>
#include <list>

namespace Verdandi
{


    /////////////////////////////////////
    // GRIDTONETWORKOBSERVATIONMANAGER //
    /////////////////////////////////////


    //! Observation operator that maps from a grid to a network.
    /*! The operator linearly interpolates from a regular grid to a list of
      locations (i.e., a network).
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class GridToNetworkObservationManager: public VerdandiBase
    {
    public:
        typedef Matrix<T, General, RowSparse> tangent_operator_matrix;
        typedef Matrix<T, General, RowSparse> error_variance;
        typedef Vector<T> tangent_operator_row;

    protected:

        /*** Observations ***/

        //! File that stores the observations.
        string observation_file_;
        //! Size in bytes of one observations vector.
        int observation_size_;
        //! Period with which observations are available.
        double Delta_t_;
        //! Period with which available observations are actually loaded.
        int Nskip_;
        //! Indicates if observation has been loaded yet.
        bool observation_loaded_;
        //! Duration during which observations are assimilated.
        double final_date_;

        //! Number total of observations date.
        int Nobservation_date_;
        //! Observation date.
        Vector<double> observation_date_;
        //! Lower bound of the last observation time interval.
        double date_inf_;
        //! Upper bound of the last observation time interval.
        double date_sup_;

        //! Number of total observations at current date.
        int Nobservation_;
        //! Observation data.
        Vector<T> observation_;

        //! Index along x.
        Vector<T> location_x_;
        //! Index along y.
        Vector<T> location_y_;

        //! Interpolation indices for all locations.
        Matrix<int> interpolation_index_;

        //! Interpolation weights for all locations.
        Matrix<T> interpolation_weight_;

        //! Interpolation indices for active locations.
        Matrix<int> active_interpolation_index_;

        //! Interpolation weights for active locations.
        Matrix<T> active_interpolation_weight_;

        //! Availability of observations at current date.
        bool availability_;

        //! Observation error variance.
        T error_variance_value_;

        /*** Model domain ***/

        //! Number of points along x in the model.
        int Nx_model_;
        //! Number of points along y in the model.
        int Ny_model_;

    public:
        // Constructors and destructor.
        GridToNetworkObservationManager();
        template <class ClassModel>
        GridToNetworkObservationManager(const ClassModel& model,
                                        string configuration_file);
        ~GridToNetworkObservationManager();

        // Initialization.
        template <class ClassModel>
        void Initialize(const ClassModel& model, string configuration_file);

        void SetAllActive();

        template <class Model>
        void SetDate(const Model& model, double date);
        template <class Model>
        void SetDate(const Model& model, double date_inf, double date_sup,
                     bool left_closed = true, bool right_closed = false);

        // Loads the observations.
        void LoadObservation();
        // Loads the observations at a given date.
        void LoadObservation(double date, Vector<T>& observation);

        bool HasObservation() const;

        // Access.
        int GetNobservation() const;
        const Vector<T>& GetObservation() const;
        bool IsOperatorSparse() const;
        bool IsErrorSparse() const;
        bool HasErrorMatrix() const;

        // Operators.
        template <class state_vector>
        void ApplyOperator(const state_vector& x, Vector<T>& y) const;

        template <class state_vector>
        void ApplyTangentOperator(const state_vector& x, Vector<T>& y) const;
        T GetTangentOperator(int i, int j) const;
        void GetTangentOperatorRow(int row, tangent_operator_row&
                                   tangent_operator_row) const;
        const tangent_operator_matrix& GetTangentOperatorMatrix()
            const;

        template <class state_vector>
        void ApplyAdjointOperator(const state_vector& x, Vector<T>& y) const;

        template <class state_vector>
        void GetInnovation(const state_vector& state,
                           Vector<T>& innovation) const;

        bool HasBLUECorrection() const;
        void GetBLUECorrection(Vector<T>& BLUE_correction) const;

        T GetObservationErrorCovariance(int i, int j) const;
        const error_variance& GetObservationErrorVariance() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_HXX
#endif
