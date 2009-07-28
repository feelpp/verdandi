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
    class GridToNetworkObservationManager
    {
    protected:

        /*** Observations ***/

        //! File that stores the observations.
        string observation_file_;
        //! Period with which observations are available.
        int period_observation_;
        //! Period with which available observations are actually loaded.
        int Nskip_;

        //! Observation data.
        Vector<T> observation_;

        //! Index along x.
        Vector<T> location_x_;
        //! Index along y.
        Vector<T> location_y_;

        //! Number of total observations at current date.
        int Nobservation_;

        //! Interpolation indices for all locations.
        Matrix<int> interpolation_index_;

        //! Interpolation weights for all locations.
        Matrix<T> interpolation_weight_;

        //! Interpolation indices for active locations.
        Matrix<int> active_interpolation_index_;

        //! Interpolation weights for active locations.
        Matrix<T> active_interpolation_weight_;

        //! Is the observation operator available in a sparse matrix?
        bool operator_sparse_;

        //! Availability of observations at current date.
        bool availability_;

        //! Observation error variance.
        T error_variance_;
        //! Is the observation error covariance matrix sparse?
        bool error_sparse_;
        //! Is the observation error covariance available in a matrix?
        bool error_matrix_availability_;

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

        // Loads the observations at a given date.
        template <class ClassModel>
        void LoadObservation(const ClassModel& model);

        bool HasObservation() const;

        // Access.
        int GetNobservation() const;
        const Vector<T>& GetObservation() const;
        bool IsOperatorSparse() const;
        bool IsErrorSparse() const;
        bool HasErrorMatrix() const;

        // Operators.
        void ApplyOperator(const Vector<T>& x, Vector<T>& y) const;

        void ApplyTangentOperator(const Vector<T>& x, Vector<T>& y) const;
        T GetTangentOperator(int i, int j) const;
        void GetTangentOperatorRow(int row,
                                   Vector<T>& tangent_operator_row) const;
        const Matrix<T, General, RowSparse>& GetTangentOperatorMatrix() const;

        void ApplyAdjointOperator(const Vector<T>& x, Vector<T>& y) const;

        void GetInnovation(const Vector<T>& state,
                           Vector<T>& innovation) const;

        bool HasBLUECorrection() const;
        void GetBLUECorrection(Vector<T>& BLUE_correction) const;

        T GetObservationErrorCovariance(int i, int j) const;
        const Matrix<T>& GetObservationErrorCovarianceMatrix() const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_GRIDTONETWORKOBSERVATIONMANAGER_HXX
#endif
