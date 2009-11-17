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


#ifndef VERDANDI_FILE_LINEAROBSERVATIONMANAGER_HXX


namespace Verdandi
{


    //////////////////////////////
    // LINEAROBSERVATIONMANAGER //
    //////////////////////////////


    //! Linear observation operator.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class LinearObservationManager: public VerdandiBase
    {
    public:
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
        typedef Matrix<T, General, RowSparse> tangent_operator_matrix;
#else
        typedef Matrix<T> tangent_operator_matrix;
#endif
        typedef Matrix<T, General, RowSparse> error_variance;
        typedef Vector<T> tangent_operator_row;

    protected:

        /*** Observations ***/

        //! File that stores the observations.
        string observation_file_;
        //! Period with which observations are available.
        int period_observation_;
        //! Period with which available observations are actually loaded.
        int Nskip_;

        //! Tangent operator matrix (H).
        tangent_operator_matrix tangent_operator_matrix_;
        //! How is defined the observation operator?
        string operator_definition_;
        //! In case of a diagonal operator.
        T operator_diagonal_value_;
        //! In case of an operator defined in a file.
        string operator_file_;

        //! Observation data.
        Vector<T> observation_;

        //! Index along x.
        Vector<T> location_x_;
        //! Index along y.
        Vector<T> location_y_;

        //! Number of total observations at current date.
        int Nobservation_;

        //! Availability of observations at current date.
        bool availability_;

        //! Observation error variance.
        T error_variance_value_;
        //! Observation error covariance matrix (R).
        error_variance error_variance_;

        /*** Model domain ***/

        //! Number of points along x in the model.
        int Nx_model_;
        //! Number of points along y in the model.
        int Ny_model_;


    public:
        // Constructors and destructor.
        LinearObservationManager();
        template <class Model>
        LinearObservationManager(const Model& model,
                                 string configuration_file);
        ~LinearObservationManager();

        // Initialization.
        template <class Model>
        void Initialize(const Model& model, string configuration_file);

        // Loads the observations at a given date.
        template <class Model>
        void LoadObservation(const Model& model);

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
        const tangent_operator_matrix& GetTangentOperatorMatrix() const;
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


#define VERDANDI_FILE_LINEAROBSERVATIONMANAGER_HXX
#endif
