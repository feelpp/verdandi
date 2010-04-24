// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton, Marc Fragu
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
#include <limits>
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
        observation_aggregator_(configuration_file)
    {
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
      \param[in] model model.
      \param[in] configuration_file configuration file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    void LinearObservationManager<T>
    ::Initialize(const Model& model, string configuration_file)
    {
        observation_aggregator_.Initialize(configuration_file);

        GetPot configuration_stream(configuration_file, "#", "\n");

        Nstate_model_ = model.GetNstate();

        configuration_stream.set_prefix("observation/");
        configuration_stream.set("File", observation_file_);
        configuration_stream.set("Type", observation_type_, "", "state");
        configuration_stream.set("Delta_t",
                                 Delta_t_,  "> 0");
        configuration_stream.set("Nskip", Nskip_, "> 0");
        configuration_stream.set("Final_date", final_date_, "",
                                 numeric_limits<double>::max());

        configuration_stream.set("error/Variance", error_variance_value_,
                                 "> 0");

        configuration_stream.set("operator/Definition",
                                 operator_definition_, "'diagonal'| 'file'");
        configuration_stream.set("operator/Diagonal_value",
                                 operator_diagonal_value_);
        configuration_stream.set("operator/File", operator_file_);

        date_ = numeric_limits<double>::min();

        /*** Building the matrices ***/


        if (operator_definition_ == "diagonal")
        {
            Nobservation_ = Nstate_model_;
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
            build_diagonal_sparse_matrix(Nstate_model_,
                                         operator_diagonal_value_,
                                         tangent_operator_matrix_);
#else
            tangent_operator_matrix_.Reallocate(Nobservation_, Nstate_model_);
            tangent_operator_matrix_.SetIdentity();
            Mlt(operator_diagonal_value_, tangent_operator_matrix_);
#endif
        }
        else if (operator_definition_ == "file")
        {
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
            throw ErrorUndefined("LinearObservationManager::Initialize()",
                                 "File definition of sparse tangent operator "
                                 "matrix is not yet supported.");
#else
            tangent_operator_matrix_.Read(operator_file_);
            if (tangent_operator_matrix_.GetN() != model.GetNstate())
                throw ErrorArgument("LinearObservationManager::Initialize()",
                                    "The number of columns of the tangent "
                                    "operator matrix ("
                                    + to_str(tangent_operator_matrix_.GetN())
                                    + ") defined in the file \"" +
                                    operator_file_
                                    + "\" is inconsistent with the"
                                    " dimension of the model state("
                                    + to_str(model.GetNstate())  + ").");
            Nobservation_ = tangent_operator_matrix_.GetM();
#endif
        }

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
                                     error_variance_);
#else
        error_variance_.Reallocate(Nobservation_, Nobservation_);
        error_variance_.SetIdentity();
        Mlt(error_variance_value_, error_variance_);
#endif

        if (observation_type_ == "state")
            Nbyte_observation_ = Nstate_model_ * sizeof(T) + sizeof(int);

        if (observation_type_ == "observation")
            Nbyte_observation_ = Nobservation_ * sizeof(T) + sizeof(int);

        int expected_file_size;
        expected_file_size = Nbyte_observation_ *
            int(final_date_ / (Delta_t_ * Nskip_));

        int file_size;
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());

#ifdef SELDON_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw IOError("LinearObservationManager"
                          "::Initialize(model, configuration_file)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif

        file_stream.seekg(0, ios_base::end);
        file_size = file_stream.tellg() ;
        file_stream.close();

        if (expected_file_size > file_size)
            throw IOError("LinearObservationManager"
                          "::Initialize(model, configuration_file)",
                          "Too few available observations, the size of \""
                          + observation_file_ + "\" must be greater than "
                          + to_str(expected_file_size) + ".");
    }


    /////////////////////////////
    // OBSERVATIONS MANAGEMENT //
    /////////////////////////////


    //! Creates a new track.
    /*!
      \return The index of the new track.
    */
    template <class T>
    int LinearObservationManager<T>
    ::CreateTrack()
    {
        return observation_aggregator_.CreateTrack();
    }


    //! Sets the track to a given track.
    /*!
      \param[in] track the given track.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetTrack(int track)
    {
        observation_aggregator_.SetTrack(track);
    }


    //! Sets the date of observations to be loaded.
    /*!
      \param[in] model the model.
      \param[in] date a given date.
    */
    template <class T>
    template <class Model>
    void LinearObservationManager<T>
    ::SetDate(const Model& model, double date)
    {
        SetDate(date);
    }


    //! Sets the date of observations to be loaded.
    /*!
      \param[in] date a given date.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetDate(double date)
    {
        if (date_ == date)
            return;

        date_ = date;
        SetAvailableDate(date_, available_date_);
    }


    //! Sets the available observation dates at a given date.
    /*!
      \param[in] date the given date.
      \param[out] available_date the available observation dates.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetAvailableDate(double date, LinearObservationManager<T>::date_vector&
                       available_date) const
    {
        double date_inf, date_sup;
        observation_aggregator_.GetContributionInterval(date, date_inf,
                                                        date_sup);
        SetAvailableDate(date_inf, date_sup, available_date);

        Logger::Log<3>(*this, to_str(date) + ", [" + to_str(date_inf) + " " +
                       to_str(date_sup) + "], {" + to_str(available_date) +
                       "}\n");
    }


    //! Sets available observation dates at a given time interval.
    /*!
      \param[in] date_inf lower bound of the given time interval.
      \param[in] date_sup upper bound of the given time interval.
      \param[out] available_date the available observation dates.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetAvailableDate(double date_inf, double date_sup,
                       LinearObservationManager<T>::date_vector&
                       available_date) const
    {
        available_date.Clear();

        date_sup = date_sup < final_date_ ? date_sup : final_date_;

        double period = Delta_t_ * Nskip_;
        double available_date_0 = floor(date_inf / period) * period;
        if (available_date_0 == date_inf)
            available_date.PushBack(available_date_0);
        available_date_0 += period;
        for (double t = available_date_0; t < date_sup; t += period)
            available_date.PushBack(t);
    }


    ////////////////////////////
    // FLATTENED OBSERVATIONS //
    ////////////////////////////


    /*** Gets observations ***/


    //! Gets observations flattened over a list of dates.
    /*! The observations available at date \a date are loaded and concatenated
      in a vector.
      \param[in] date the given date.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetFlattenedObservation(available_date, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup]
      are loaded and concatenated in a vector.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date_inf, double date_sup,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetFlattenedObservation(available_date, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(
        LinearObservationManager<T>::observation_vector& observation)
    {
        GetFlattenedObservation(available_date_, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and concatenated in a vector.
      \param[in] available_date the given observation date vector.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(const LinearObservationManager<T>::date_vector&
                              available_date,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        observation_vector2 observation2;
        GetRawObservation(available_date, observation2);
        observation2.Flatten(observation);
    }


    /*** Gets observations and associated variables ***/


    //! Gets observations flattened over a list of dates.
    /*! The observations available at date \a date are loaded and concatenated
      in a vector.
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date,
                              LinearObservationManager<T>::variable_vector&
                              observation_variable,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetFlattenedObservation(available_date, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup]
      are loaded and concatenated in a vector.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date_inf, double date_sup,
                              LinearObservationManager<T>::variable_vector&
                              observation_variable,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetFlattenedObservation(available_date, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::observation_vector& observation)
    {
        GetFlattenedObservation(available_date_, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and concatenated in a vector.
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(const LinearObservationManager<T>::date_vector&
                              available_date,
                              LinearObservationManager<T>::variable_vector&
                              observation_variable,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        GetRawObservation(available_date, observation_variable2,
                          observation3);
        observation3.Flatten(observation);
        observation_variable2.Flatten(observation_variable);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets observations flattened over a list of dates.
    /*! The observations available at date \a date are loaded and concatenated
      in a vector.
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date,
                              LinearObservationManager<T>::variable_vector&
                              observation_variable,
                              LinearObservationManager<T>::index_vector&
                              observation_index,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetFlattenedObservation(available_date, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup]
      are loaded and concatenated in a vector.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double date_inf, double date_sup,
                              LinearObservationManager<T>::variable_vector&
                              observation_variable,
                              LinearObservationManager<T>::index_vector&
                              observation_index,
                              LinearObservationManager<T>::observation_vector&
                              observation)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetFlattenedObservation(available_date, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector& observation_index,
        LinearObservationManager<T>::observation_vector& observation)
    {
        GetFlattenedObservation(available_date_, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and concatenated in a vector.
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector& observation_index,
        LinearObservationManager<T>::observation_vector& observation)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        index_vector3 observation_index3;
        GetRawObservation(available_date, observation_variable2,
                          observation_index3, observation3);
        observation3.Flatten(observation);
        observation_variable2.Flatten(observation_variable);
        observation_index3.Flatten(observation_index);
    }


    /////////////////////////////
    // AGGREGATED OBSERVATIONS //
    /////////////////////////////


    /*** Gets observations ***/


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the given date are loaded and
      aggregated.
      \param[in] date the given date.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date,
        LinearObservationManager<T>::observation_vector& observation)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetAggregatedObservation(available_date, observation);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup] are loaded
      and aggregated.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::observation_vector& observation)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetAggregatedObservation(available_date, observation);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        LinearObservationManager<T>::observation_vector& observation)
    {
        GetAggregatedObservation(available_date_, observation);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and aggregated.
      \param[in] available_date the given observation date vector.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::observation_vector& observation)
    {
        observation_vector2 observation2;
        GetRawObservation(available_date, observation2);
        observation_aggregator_.Aggregate(available_date, observation2, date_,
                                          observation);
    }


    /*** Gets observations and associated variables ***/


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the given date are loaded and
      aggregated.
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetAggregatedObservation(available_date, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup] are loaded
      and aggregated.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetAggregatedObservation(available_date, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        GetAggregatedObservation(available_date_, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and aggregated.
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        GetRawObservation(available_date, observation_variable2,
                          observation3);
        observation_aggregator_.Aggregate(available_date,
                                          observation_variable2,
                                          observation3,
                                          date_,
                                          observation_variable,
                                          observation2);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the given date are loaded and
      aggregated.
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector2& observation_index2,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetAggregatedObservation(available_date, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations in the interval [\a date_inf, \a date_sup] are loaded
      and aggregated.
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector2& observation_index2,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetAggregatedObservation(available_date, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector2& observation_index2,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        GetAggregatedObservation(available_date_, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of dates.
    /*! The observations available at the dates \a available_date are loaded
      and aggregated.
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector& observation_variable,
        LinearObservationManager<T>::index_vector2& observation_index2,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        index_vector3 observation_index3;
        GetRawObservation(available_date, observation_variable2,
                          observation_index3, observation3);
        observation_aggregator_.Aggregate(available_date,
                                          observation_variable2,
                                          observation_index3,
                                          observation3,
                                          date_,
                                          observation_variable,
                                          observation_index2,
                                          observation2);
    }


    //////////////////
    // OBSERVATIONS //
    //////////////////


    /*** Gets observations ***/


    //! Gets available observations at a given date.
    /*!
      \param[in] date the given date.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetRawObservation(available_date, observation2);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetRawObservation(available_date, observation2);
    }


    //! Gets observations at the current date.
    /*!
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        GetRawObservation(available_date_, observation2);
    }


    //! Gets observations of a list of dates.
    /*!
      \param[in] available_date the given observation date vector.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::observation_vector2& observation2)
    {
        ReadObservation(available_date, observation2);
    }


    /*** Gets observations and associated variables ***/


    //! Gets available observations at a given date.
    /*!
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetRawObservation(available_date, observation_variable2,
                          observation3);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetRawObservation(available_date, observation_variable2,
                          observation3);
    }



    //! Gets observations at the current date.
    /*!
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        GetRawObservation(available_date_, observation_variable2,
                          observation3);
    }


    //! Gets observations of a list of dates.
    /*!
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        ReadObservationVariable(available_date, observation_variable2);
        ReadObservation(available_date, observation_variable2, observation3);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets available observations at a given date.
    /*!
      \param[in] date the given date.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::index_vector3& observation_index3,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        date_vector& available_date;
        SetAvailableDate(date, available_date);
        GetRawObservation(available_date, observation_variable2,
                          observation_index3, observation3);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] date_inf lower_bound of the given interval.
      \param[in] date_sup upper_bound of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double date_inf, double date_sup,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::index_vector3& observation_index3,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        date_vector& available_date;
        SetAvailableDate(date_inf, date_sup, available_date);
        GetRawObservation(available_date, observation_variable2,
                          observation_index3, observation3);
    }



    //! Gets observations at the current date.
    /*!
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::index_vector3& observation_index3,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        GetRawObservation(available_date_, observation_variable2,
                          observation_index3, observation3);
    }


    //! Gets observations of a list of dates.
    /*!
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector2& observation_variable2,
        LinearObservationManager<T>::index_vector3& observation_index3,
        LinearObservationManager<T>::observation_vector3& observation3)
    {
        ReadObservationVariable(available_date, observation_variable2);
        ReadObservation(available_date, observation_variable2, observation3);
        ReadObservationIndex(available_date, observation_variable2,
                             observation_index3);
    }


    ///////////////////////////
    // READ OBSERVATION DATA //
    ///////////////////////////


    //! Builds variables vector associated with given observations.
    /*!
      \param[in] available_date the given observation date vector.
      \param[out] observation_variable variables associated with the
      observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservationVariable(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::variable_vector2& observation_variable2)
        const
    {
        int Nt = available_date.GetSize();
        observation_variable2.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
            observation_variable2(h).PushBack(0);
    }


    //! Builds observations associated with given dates and variables.
    /*!
      \param[in] available_date the given observation date vector.
      \param[in] observation_variable variables associated with the
      observations.
      \param[out] observation the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        const LinearObservationManager<T>::variable_vector2&
        observation_variable2,
        LinearObservationManager<T>::observation_vector3& observation3) const
    {
        int Nvariable, Nt;
        Nt = available_date.GetSize();
        observation3.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
        {
            Nvariable = observation_variable2(h).GetSize();
            observation3.Reallocate(h, Nvariable);
            for (int v = 0; v < Nvariable; v++)
                ReadObservation(available_date(h),
                                observation_variable2(h, v),
                                observation3.GetVector(h, v));
        }
    }


    //! Builds observations associated with given dates.
    /*!
      \param[in] available_date the given observation date vector.
      \param[out] observation the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(
        const LinearObservationManager<T>::date_vector& available_date,
        LinearObservationManager<T>::observation_vector2& observation2) const
    {
        int Nt = available_date.GetSize();
        observation2.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
            ReadObservation(available_date(h), 0, observation2(h));
    }


    //! Reads observation from observation file given a date and a variable.
    /*!
      \param[in] date the date.
      \param[in] variable the variable.
      \param[out] observation the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(
        double date, int variable,
        LinearObservationManager<T>::observation_vector& observation) const
    {
        observation.Reallocate(Nobservation_);
        observation_vector input_data;
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());
#ifdef SELDON_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw IOError("LinearObservationManager"
                          "::LoadObservation(model)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif
        streampos position;
        position =  (int(date / (Delta_t_ * Nskip_)) + variable)
            * Nbyte_observation_;
        file_stream.seekg(position);
        input_data.Read(file_stream);
        file_stream.close();

        if (observation_type_ == "state")
            ApplyOperator(input_data, observation);
    }


    //! Reads observations indexes.
    /*!
      \param[in] available_date the available date.
      \param[in] observation_variable variable associated with the
      observations.
      \param[out] observation_index the indexes associated with the
      observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservationIndex(
        const LinearObservationManager<T>::date_vector& available_date,
        const LinearObservationManager<T>::variable_vector2&
        observation_variable2,
        LinearObservationManager<T>::index_vector3& observation_index3) const
    {
        int Nvariable, Nt;
        Nt = available_date.GetSize();
        observation_index3.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
        {
            Nvariable = observation_variable2.GetSize(h);
            observation_index3.Reallocate(h, Nvariable);
            for (int v = 0; v < Nvariable; v++)
            {
                observation_index3.Reallocate(h, v, Nobservation_);
                observation_index3(h, v).Fill(0);
            }
        }
    }


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Gets observation.
    /*!
      \param[out] observation observation vector.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetObservation(
        LinearObservationManager<T>::observation_vector& observation)
    {
        GetAggregatedObservation(observation);
    }


    ////////////////
    // INNOVATION //
    ////////////////


    //! Gets innovation.
    /*!
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class T>
    template <class state_vector>
    void LinearObservationManager<T>
    ::GetInnovation(
        const state_vector& state,
        LinearObservationManager<T>::observation_vector& innovation)
    {
        innovation.Reallocate(Nobservation_);
        observation_vector observation;
        GetObservation(observation);
        ApplyOperator(state, innovation);
        Mlt(T(-1), innovation);
        Add(T(1), observation, innovation);
    }


    ////////////
    // ACCESS //
    ////////////


    //! Indicates if some observations are available at a given date.
    /*!
      \param[in] date a given date.
    */
    template <class T>
    bool LinearObservationManager<T>::HasObservation(double date) const
    {
        if (date == date_)
            return HasObservation();
        date_vector available_date;
        SetDate(date, available_date);
        return available_date.GetSize() != 0;
    }


    //! Indicates if some observations are available at current time.
    template <class T>
    bool LinearObservationManager<T>::HasObservation() const
    {
        return available_date_.GetSize() != 0;
    }


    //! Gets Nobservation_ value.
    /*!
      \return The total number of observation at current date.
    */
    template <class T>
    int LinearObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
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
    ::ApplyOperator(const state_vector& x,
                    LinearObservationManager<T>::observation_vector& y) const
    {
        if (operator_definition_ == "diagonal")
        {
            Copy(x, y);
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
    ::ApplyTangentOperator(const state_vector& x,
                           LinearObservationManager<T>::observation_vector& y)
        const
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
        if (operator_definition_ == "diagonal")
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
        if (operator_definition_ == "diagonal")
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
        return tangent_operator_matrix_;
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    template <class state_vector>
    void LinearObservationManager<T>
    ::ApplyAdjointOperator(
        const state_vector& x,
        LinearObservationManager<T>::observation_vector& y)
        const
    {
        if (operator_definition_ == "diagonal")
        {
            y = x;
            Mlt(operator_diagonal_value_, y);
        }

        // Operator defined in a file.
        else
            Mlt(1., SeldonTrans, tangent_operator_matrix_, x, 0., y);
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
        return error_variance_;
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
    }


} // namespace Verdandi.


#define VERDANDI_FILE_LINEAROBSERVATIONMANAGER_CXX
#endif
