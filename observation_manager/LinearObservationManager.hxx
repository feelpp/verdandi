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

#include <iostream>
#include <list>

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

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        typedef Matrix<T, General, RowSparse> error_variance;
#else
        typedef Matrix<T> error_variance;
#endif
        typedef Vector<T> tangent_operator_row;

    protected:

        /*** Observations ***/

        //! File that stores the observations.
        string observation_file_;
        //! How are stored the observations.
        string observation_type_;
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

        //! Observation data.
        Vector<T> observation_;
        //! Observation date.
        Vector<double> observation_date_;
        //! Lower bound of the last observation time interval.
        double date_inf_;
        //! Upper bound of the last observation time interval.
        double date_sup_;

        //! Number total of observations at current date.
        int Nobservation_;
        //! Number total of observations date.
        int Nobservation_date_;

        //! Availability of observations at current date.
        bool availability_;

        //! Tangent operator matrix (H).
        tangent_operator_matrix tangent_operator_matrix_;
        //! How is defined the observation operator?
        string operator_definition_;
        //! In case of a diagonal operator.
        T operator_diagonal_value_;
        //! In case of an operator defined in a file.
        string operator_file_;


        //! Observation error variance.
        T error_variance_value_;
        //! Observation error covariance matrix (R).
        error_variance error_variance_;

        /*** Model domain ***/

        int Nstate_model_;


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
        void SetDate(const Model& model, double date);
        // Loads the observations at a given date.
        template <class Model>
        void SetDate(const Model& model, double date_inf, double date_sup,
                     bool left_closed = true, bool right_closed = false);
        // Loads the observations at a given date.
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
