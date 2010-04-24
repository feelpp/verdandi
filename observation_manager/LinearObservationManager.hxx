// Copyright (C) 2008-2009 INRIA
// Author(s): Claire Mouton, Vivien Mallet, Marc Fragu
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

#ifndef OBSERVATION_AGGREGATOR
#define OBSERVATION_AGGREGATOR ObservationAggregator
#endif

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#include QUOTE(OBSERVATION_AGGREGATOR.cxx)


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
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, RowSparse> tangent_operator_matrix;
#else
        //! Type of the tangent linear operator.
        typedef Matrix<T> tangent_operator_matrix;
#endif

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        //! Type of the observation error covariance matrix.
        typedef Matrix<T, General, RowSparse> error_variance;
#else
        //! Type of the observation error covariance matrix.
        typedef Matrix<T> error_variance;
#endif
        //! Type of a row of the tangent linear operator.
        typedef Vector<T> tangent_operator_row;

        typedef Vector<T> observation_vector;
        typedef Vector2<T> observation_vector2;
        typedef Vector3<T> observation_vector3;

        typedef Vector<int> variable_vector;
        typedef Vector2<int> variable_vector2;
        typedef Vector3<int> variable_vector3;

        typedef Vector<int> index_vector;
        typedef Vector2<int> index_vector2;
        typedef Vector3<int> index_vector3;

        typedef Vector<double> date_vector;
        typedef Vector2<double> date_vector2;
        typedef Vector3<double> date_vector3;

    protected:

        /*** Observation file structure ***/

        //! File that stores the observations.
        string observation_file_;
        //! How are stored the observations.
        string observation_type_;
        //! Total number of observations at current date.
        int Nobservation_;
        //! Size in bytes of an observation vector.
        size_t Nbyte_observation_;
        //! Period with which observations are available.
        double Delta_t_;
        //! Period with which available observations are actually loaded.
        int Nskip_;
        //! Duration during which observations are assimilated.
        double final_date_;

        /*** Observation dates ***/

        //! Requested date.
        double date_;
        //! Available observation date of the time interval.
        date_vector available_date_;
        //! Observations aggregator.
        OBSERVATION_AGGREGATOR<T> observation_aggregator_;

        /*** Observation operator ***/

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

        //! The size of a model state.
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

        int CreateTrack();
        void SetTrack(int track);

        template <class Model>
        void SetDate(const Model& model, double date);
        void SetDate(double date);
        void SetAvailableDate(double date, date_vector& available_date) const;
        void SetAvailableDate(double date_inf, double date_sup,
                              date_vector& available_date) const;


        ////////////////////////////
        // FLATTENED OBSERVATIONS //
        ////////////////////////////


        /*** Gets observations ***/

        void GetFlattenedObservation(double date,
                                     observation_vector& observation);
        void GetFlattenedObservation(double date_inf, double date_sup,
                                     observation_vector& observation);
        void GetFlattenedObservation(observation_vector& observation);
        void GetFlattenedObservation(const date_vector& available_date,
                                     observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetFlattenedObservation(double date,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(double date_inf, double date_sup,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(const date_vector& available_date,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);

        /*** Gets observations, associated variables and associated index ***/

        void GetFlattenedObservation(double date,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(double date_inf, double date_sup,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(const date_vector& available_date,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);


        /////////////////////////////
        // AGGREGATED OBSERVATIONS //
        /////////////////////////////


        /*** Gets observations ***/

        void GetAggregatedObservation(double date,
                                      observation_vector& observation);
        void GetAggregatedObservation(double date_inf, double date_sup,
                                      observation_vector& observation);
        void GetAggregatedObservation(observation_vector& observation);
        void GetAggregatedObservation(const date_vector& available_date,
                                      observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetAggregatedObservation(double date,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double date_inf, double date_sup,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const date_vector& available_date,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);

        /*** Gets observations, associated variables and associated index ***/

        void GetAggregatedObservation(double date,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double date_inf, double date_sup,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const date_vector& available_date,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);


        //////////////////////
        // RAW OBSERVATIONS //
        //////////////////////


        /*** Gets observations ***/

        void GetRawObservation(double date,
                               observation_vector2& observation2);
        void GetRawObservation(double date_inf, double date_sup,
                               observation_vector2& observation2);
        void GetRawObservation(observation_vector2& observation2);
        void GetRawObservation(const date_vector& available_date,
                               observation_vector2& observation2);

        /*** Gets observations and associated variables ***/

        void GetRawObservation(double date,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(double date_inf, double date_sup,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(const date_vector& available_date,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);

        /*** Gets observations, associated variables and associated index ***/

        void GetRawObservation(double date,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(double date_inf, double date_sup,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(const date_vector& available_date,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);


        ///////////////////////////////
        // READ OBSERVATIONS METHODS //
        ///////////////////////////////


        void ReadObservationVariable(const date_vector& available_date,
                                     variable_vector2& observation_variable2)
            const;
        void ReadObservation(const date_vector& available_date,
                             const variable_vector2& observation_variable2,
                             observation_vector3& observation3) const;
        void ReadObservation(const date_vector& available_date,
                             observation_vector2& observation2) const;
        void ReadObservation(double date, int variable,
                             observation_vector& observation) const;
        void ReadObservationIndex(const date_vector& available_date, const
                                  variable_vector2& observation_variable2,
                                  index_vector3& observation_index3) const;


        /////////////////
        // OBSERVATION //
        /////////////////


        void GetObservation(observation_vector& observation);


        ////////////////
        // INNOVATION //
        ////////////////


        template <class state_vector>
        void GetInnovation(const state_vector& state,
                           observation_vector& innovation);


        ////////////
        // ACCESS //
        ////////////


        bool HasObservation() const;
        bool HasObservation(double date) const;
        int GetNobservation() const;
        bool IsOperatorSparse() const;
        bool IsErrorSparse() const;
        bool HasErrorMatrix() const;


        ///////////////
        // OPERATORS //
        ///////////////


        template <class state_vector>
        void ApplyOperator(const state_vector& x, observation_vector& y)
            const;

        template <class state_vector>
        void ApplyTangentOperator(const state_vector& x,
                                  observation_vector& y) const;
        T GetTangentOperator(int i, int j) const;
        void GetTangentOperatorRow(int row, tangent_operator_row&
                                   tangent_operator_row) const;
        const tangent_operator_matrix& GetTangentOperatorMatrix() const;
        template <class state_vector>
        void ApplyAdjointOperator(const state_vector& x,
                                  observation_vector& y) const;

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
