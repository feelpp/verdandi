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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_LINEAROBSERVATIONMANAGER_HXX

#include <iostream>
#include <list>

#ifndef OBSERVATION_AGGREGATOR
#define OBSERVATION_AGGREGATOR ObservationAggregator
#endif

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#include QUOTE(OBSERVATION_AGGREGATOR.hxx)


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

#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
#ifdef VERDANDI_WITH_MPI
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, PETScMPIAIJ> tangent_linear_operator;
#else
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, RowSparse> tangent_linear_operator;
#endif
#else
#ifdef VERDANDI_WITH_MPI
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, PETScMPIDense> tangent_linear_operator;
#else
        //! Type of the tangent linear operator.
        typedef Matrix<T> tangent_linear_operator;
#endif
#endif


#ifdef VERDANDI_ERROR_SPARSE
#ifdef VERDANDI_WITH_MPI
        //! Type of the observation error covariance matrix.
        typedef Matrix<T, General, PETScMPIAIJ> error_variance;
#else
        //! Type of the observation error covariance matrix.
        typedef Matrix<T, General, RowSparse> error_variance;
#endif
#else
#ifdef VERDANDI_WITH_MPI
        //! Type of the observation error covariance matrix.
        typedef Matrix<T, General, PETScMPIDense> error_variance;
#else
        //! Type of the observation error covariance matrix.
        typedef Matrix<T> error_variance;
#endif
#endif

#ifdef VERDANDI_WITH_MPI
        //! Type of a row of the tangent linear operator.
        typedef Vector<T, PETScPar> tangent_linear_operator_row;

        //! Type of the observation vector.
        typedef Vector<T, PETScPar> observation;
        //! Type of the observation vector.
        typedef Vector<T, PETScPar> observation_vector;

#else
        //! Type of a row of the tangent linear operator.
        typedef Vector<T> tangent_linear_operator_row;

        //! Type of the observation vector.
        typedef Vector<T> observation;
        //! Type of the observation vector.
        typedef Vector<T> observation_vector;
#endif

        //! Type of the observation vector 2.
        typedef Vector2<T> observation_vector2;
        //! Type of the observation vector 3.
        typedef Vector3<T> observation_vector3;

        //! Type of the variable vector.
        typedef Vector<int> variable_vector;
        //! Type of the variable vector 2.
        typedef Vector2<int> variable_vector2;
        //! Type of the variable vector 3.
        typedef Vector3<int> variable_vector3;

        //! Type of the index vector.
        typedef Vector<int> index_vector;
        //! Type of the index vector 2.
        typedef Vector2<int> index_vector2;
        //! Type of the index vector 3.
        typedef Vector3<int> index_vector3;

        //! Type of the time vector.
        typedef Vector<double> time_vector;
        //! Type of the time vector 2.
        typedef Vector2<double> time_vector2;
        //! Type of the time vector 3.
        typedef Vector3<double> time_vector3;

    protected:

        /*** Observation file structure ***/

        //! File that stores the observations.
        string observation_file_;
        //! Type of the file.
        string observation_file_type_;
        //! Path to the dataset where observations are stored (HDF5 filetype).
        string observation_dataset_path_;
        //! How are stored the observations.
        string observation_type_;
        //! Total number of observations at current time.
        size_t Nobservation_;
        //! Size in bytes of an observation vector.
        size_t Nbyte_observation_;
        //! Is the time interval between two observations constant?
        bool is_delta_t_constant_;
        //! Period with which observations are available (if constant).
        double Delta_t_;
        /*! Times at which observations are available (if the time interval
          between two observations is not constant). */
        Vector<double> observation_time_;

        /*! Path to the observation times (needed if the time interval between
          two observations is not constant). */
        string observation_time_file_;
        //! Period with which available observations are actually loaded.
        int Nskip_;
        //! First time at which observations are available.
        double initial_time_;
        //! Final time at which observations are available.
        double final_time_;

        /*** Observation times ***/

        //! Requested time.
        double time_;
        //! Available observation time of the time interval.
        time_vector available_time_;
        //! Contribution associated with available observations.
        Vector<double> contribution_;

        //! Observations aggregator.
        OBSERVATION_AGGREGATOR<T> observation_aggregator_;

        /*** Observation operator ***/

        //! Tangent operator matrix (H).
        tangent_linear_operator tangent_operator_matrix_;
        //! Is the operator a scaled identity matrix?
        bool operator_scaled_identity_;
        //! In case of a scaled identity operator.
        T operator_diagonal_value_;

        //! Observation error variance.
        T error_variance_value_;
        //! Observation error covariance matrix (R).
        error_variance error_variance_;
        //! Inverse of the observation error covariance matrix (R).
        error_variance error_variance_inverse_;

        /*** Triangle interpolation ***/

        //! File that stores the observations.
        string width_file_;

        /*** Model domain ***/

        //! The size of a model state.
        size_t Nstate_model_;

        //! Observation currently stored.
        observation observation_;
        //! Innovation currently stored.
        observation innovation_;
        //! Index of the row of H currently stored.
        int current_row_;
        //! Value of the row of H currently stored.
        tangent_linear_operator_row tangent_operator_row_;

    public:
        // Constructors and destructor.
        LinearObservationManager();
        template <class Model>
        LinearObservationManager(Model& model,
                                 string configuration_file);
        ~LinearObservationManager();

        // Initialization.
        template <class Model>
        void Initialize(Model& model, string configuration_file);
        template <class Model>
        void InitializeOperator(Model& model,
                                string configuration_file);

        void DiscardObservation(bool discard_observation);
        int CreateTrack();
        void SetTrack(int track);

        template <class Model>
        void SetTime(Model& model, double time);
        void SetTime(double time);
        void SetAvailableTime(double time, time_vector& available_time);
        void SetAvailableTime(double time_inf, double time_sup, time_vector&
                              available_time);
        void SetAvailableTime(double time, double time_inf, double time_sup,
                              int selection_policy,
                              time_vector& available_time);


        ////////////////////////////
        // FLATTENED OBSERVATIONS //
        ////////////////////////////


        /*** Gets observations ***/

        void GetFlattenedObservation(double time,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     observation_vector& observation);
        void GetFlattenedObservation(observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetFlattenedObservation(double time,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);

        /*** Gets observations, associated variables and associated index ***/

        void GetFlattenedObservation(double time,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);


        /////////////////////////////
        // AGGREGATED OBSERVATIONS //
        /////////////////////////////


        /*** Gets observations ***/

        void GetAggregatedObservation(double time,
                                      observation_vector& observation);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      observation_vector& observation);
        void GetAggregatedObservation(observation_vector& observation);
        void GetAggregatedObservation(const time_vector& available_time,
                                      observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetAggregatedObservation(double time,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const time_vector& available_time,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);

        /*** Gets observations, associated variables and associated index ***/

        void GetAggregatedObservation(double time,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const time_vector& available_time,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);


        //////////////////////
        // RAW OBSERVATIONS //
        //////////////////////


        /*** Gets observations ***/

        void GetRawObservation(double time,
                               observation_vector2& observation2);
        void GetRawObservation(double time_inf, double time_sup,
                               observation_vector2& observation2);
        void GetRawObservation(observation_vector2& observation2);
        void GetRawObservation(const time_vector& available_time,
                               observation_vector2& observation2);

        /*** Gets observations and associated variables ***/

        void GetRawObservation(double time,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(double time_inf, double time_sup,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(const time_vector& available_time,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);

        /*** Gets observations, associated variables and associated index ***/

        void GetRawObservation(double time,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(double time_inf, double time_sup,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(const time_vector& available_time,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);


        ///////////////////////////////
        // READ OBSERVATIONS METHODS //
        ///////////////////////////////


        void ReadObservationVariable(const time_vector& available_time,
                                     variable_vector2& observation_variable2)
            const;
        void ReadObservation(const time_vector& available_time,
                             const variable_vector2& observation_variable2,
                             observation_vector3& observation3) const;
        void ReadObservation(const time_vector& available_time,
                             observation_vector2& observation2) const;
        void ReadObservation(ifstream& file_stream, double time, int variable,
                             observation_vector& observation) const;
#ifdef VERDANDI_WITH_HDF5
        void ReadObservation(ifstream& file_stream, double time, int variable,
                             string dataset_path,
                             observation_vector& observation) const;
#endif
        void ReadObservationIndex(const time_vector& available_time, const
                                  variable_vector2& observation_variable2,
                                  index_vector3& observation_index3) const;
        void ReadObservationTriangleWidth(double time_inf, double time_sup,
                                          Vector<double>& width_left,
                                          Vector<double>& width_right) const;


        /////////////////
        // OBSERVATION //
        /////////////////


        observation& GetObservation();


        ////////////////
        // INNOVATION //
        ////////////////


        template <class state>
        observation& GetInnovation(const state& x);


        ////////////
        // ACCESS //
        ////////////


        bool HasObservation() const;
        int GetNobservation() const;
        bool IsOperatorSparse() const;
        bool IsErrorSparse() const;
        bool HasErrorMatrix() const;
        template <class state, class mat>
        void GetNudgingMatrix(const state& x, mat& M) const;

        ///////////////
        // OPERATORS //
        ///////////////


        template <class state>
        void ApplyOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyTangentLinearOperator(const state& x, observation& y) const;
        T GetTangentLinearOperator(int i, int j) const;
        tangent_linear_operator_row& GetTangentLinearOperatorRow(int row);
        const tangent_linear_operator& GetTangentLinearOperator() const;
        template <class state>
        void ApplyAdjointOperator(const state& x, observation& y) const;

        bool HasBLUECorrection() const;
        void GetBLUECorrection(observation& BLUE_correction) const;

        T GetErrorVariance(int i, int j) const;
        const error_variance& GetErrorVariance() const;
        const error_variance& GetErrorVarianceInverse() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_LINEAROBSERVATIONMANAGER_HXX
#endif
