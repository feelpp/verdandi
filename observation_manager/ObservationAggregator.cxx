// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_OSERVATIONAGGREGATOR_CXX


#include "ObservationAggregator.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    template <class T>
    ObservationAggregator<T>::ObservationAggregator()
    {
    }


    //! Constructor.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    ObservationAggregator<T>::ObservationAggregator(string configuration_file)
    {
    }


    //! Destructor.
    template <class T>
    ObservationAggregator<T>::~ObservationAggregator()
    {
    }


    //! Initializer.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ObservationAggregator<T>::Initialize(string configuration_file)
    {
        GetPot configuration_stream(configuration_file);

        configuration_stream.
            set_prefix("observation/observation_aggregator/");

        string interpolation_type_str;
        configuration_stream.set("Type", interpolation_type_str, "", "Step");
        if (interpolation_type_str == "Step")
            interpolation_type_ = type_step_;
        else if (interpolation_type_str == "Triangle")
            interpolation_type_ = type_triangle_;

        configuration_stream.set("Width_left", width_left_);
        configuration_stream.set("Width_right", width_right_);

        configuration_stream.set("Discard_observation", discard_observation_);

        active_track_index_ = CreateTrack();
    }


    //! Returns the contribution time interval corresponding to a given date.
    /*! This method returns the time interval into which observations have a
      non-zero contribution at date \a date.
      \param[in] date a given date.
      \param[out] date_inf lower bound of the time interval.
      \param[out] date_sup upper bound of the time interval.
    */
    template <class T>
    void ObservationAggregator<T>
    ::GetContributionInterval(double date, double& date_inf, double& date_sup)
        const
    {
        //! Always switch off observation.
        if (discard_observation_)
        {
            if (date_.GetSize() == 0)
            {
                date_inf = date - width_left_;
                date_inf = date_inf > 0. ? date_inf : 0.;
                date_sup = date + width_right_;
                return;
            }

            FindInterval(date_(active_track_index_), date, date_inf,
                         date_sup);

            if (date_inf == -1.)
            {
                date_inf = date - width_left_;
                date_inf = date_inf > 0. ? date_inf : 0.;
                date_sup = date + width_right_ < date_(active_track_index_)(0)
                    - width_left_ ? date + width_right_ :
                    date_(active_track_index_)(0) - width_left_;
                return;
            }

            if (date_sup == -1.)
            {
                date_inf = date - width_left_ > date_inf + width_right_ ?
                    date - width_left_ : date_inf + width_right_;
                date_inf = date_inf > 0. ? date_inf : 0.;
                date_sup = date + width_right_;
                return;
            }

            date_inf = date - width_left_ > date_inf + width_right_ ?
                date - width_left_ : date_inf + width_right_;
            date_inf = date_inf > 0. ? date_inf : 0.;
            date_sup = date + width_right_ < date_sup - width_left_ ?
                date + width_right_ : date_sup - width_left_;
            return;
        }

        date_inf = date - width_left_;
        date_inf = date_inf > 0. ? date_inf : 0.;
        date_sup = date + width_right_;
    }


    //! Computes an aggregated observation vector over a list of dates.
    /*! The observations that have a non-zero contribution at date \a date are
      aggregated.
      \param[in] observation_date the dates of \a observation.
      \param[in] observation the observations to be aggregated.
      \param[in] date the date at which the observations should be aggregated.
      \param[out] aggregated_observation the aggregated observation vector.
    */
    template <class T>
    template <class date_vector, class observation_vector2,
              class observation_vector>
    void ObservationAggregator<T>
    ::Aggregate(const date_vector& observation_date,
                const observation_vector2& observation,
                double date,
                observation_vector& aggregated_observation)
    {

        /*** Computes contribution ***/

        Vector<double> contribution;
        Contribution(date, observation_date, contribution);

        double sum = Norm1(contribution);
        if (sum == 0.)
            throw ErrorProcessing("ObservationAggregator::Aggregate()",
                                  "The sum of contributions should not be"
                                  " zero.");

        /*** Computes aggregated observations ***/

        // Assumes 'Nobservation' is constant.
        int Nobservation = observation(0).GetSize();
        aggregated_observation.Reallocate(Nobservation);
        aggregated_observation.Fill(T(0.));
        int Nt = observation_date.GetSize();
        for (int h = 0; h < Nt; h++)
            Add(contribution(h), observation(h), aggregated_observation);
        Mlt(T(1. / sum), aggregated_observation);

        if (discard_observation_)
            PushDate(date);
    }


    //! Computes aggregated observation over a list of dates and observations.
    /*! The observations that have a non-zero contribution at date \a date are
      aggregated. The variables associated with the new aggregated
      observations vector are stored in \a aggregated_variable.
      \param[in] observation_date the date of the given observations.
      \param[in] observation_variable variables associated with the
      observations.
      \param[in] observation the given observation observation.
      \param[in] date a given date.
      \param[out] aggregated_variable the variables associated with the
      aggregated observations.
      \param[out] aggregated_observation the aggregated observation.
    */
    template <class T>
    template <class date_vector, class variable_vector2,
              class observation_vector3,
              class variable_vector, class observation_vector2>
    void ObservationAggregator<T>
    ::Aggregate(const date_vector& observation_date,
                const variable_vector2& observation_variable,
                const observation_vector3& observation,
                double date,
                variable_vector& aggregated_variable,
                observation_vector2& aggregated_observation)
    {
        // Assumes 'Nobservation' is constant.
        int Nobservation = observation(0, 0).GetSize();

        /*** Builds variable vector ***/

        observation_variable.Flatten(aggregated_variable);
        RemoveDuplicate(aggregated_variable);
        int Nvariable = aggregated_variable.GetSize();

        map<int, int> variable_index_map;
        for (int v = 0; v < Nvariable; v++)
            variable_index_map[aggregated_variable(v)] = v;

        /*** Computes contribution ***/

        Vector<double> contribution;
        Contribution(date, observation_date, contribution);

        Vector<double> sum(Nvariable);
        sum.Fill(0.);

        /*** Builds aggregated observations ***/

        aggregated_observation.Reallocate(Nvariable);
        for (int v = 0; v < Nvariable; v++)
            aggregated_observation(v).Reallocate(Nobservation);
        aggregated_observation.Fill(T(0.));

        int Nt, variable_index;
        Nt = observation_date.GetSize();
        for (int h = 0; h < Nt; h++)
            for (int v = 0; v < observation(h).GetSize(); v++)
            {
                variable_index =
                    variable_index_map[observation_variable(h, v)];
                sum(variable_index) += contribution(h);
                Add(contribution(h), observation(h, v),
                    aggregated_observation(variable_index));
            }

        for (int v = 0; v < Nvariable; v++)
        {
            if (sum(v) == 0.)
                throw ErrorProcessing("ObservationAggregator::Aggregate()",
                                      "The sum of contributions for variable "
                                      + to_str(v) + " should not be zero.");
            Mlt(T(1. / sum(v)), aggregated_observation(v));
        }

        if (discard_observation_)
            PushDate(date);
    }


    //! Computes aggregated observation over a list of dates and observations.
    /*! The observations that have a non-zero contribution at date \a date are
      aggregated. The variables associated with the new aggregated
      observations vector are stored in \a aggregated_variable. The indexes
      associated with the new aggregated observations vector are stored in
      \a aggregated_index.
      \param[in] observation_date the date of the given observations.
      \param[in] observation_variable variables associated with observations.
      \param[in] observation_index corresponding observation locations.
      \param[in] observation the given observation observation.
      \param[in] date a given date.
      \param[out] aggregated_variable variables associated with the
      aggregated observations.
      \param[out] aggregated_index the aggregated locations.
      \param[out] aggregated_observation the aggregated observation.
    */
    template <class T>
    template <class date_vector, class variable_vector2,
              class index_vector3, class observation_vector3,
              class variable_vector, class index_vector2,
              class observation_vector2>
    void ObservationAggregator<T>
    ::Aggregate(const date_vector& observation_date,
                const variable_vector2& observation_variable,
                const index_vector3& observation_index,
                const observation_vector3& observation,
                double date,
                variable_vector& aggregated_variable,
                index_vector2& aggregated_index,
                observation_vector2& aggregated_observation)
    {
        throw ErrorUndefined("ObservationAggregator::Aggregate"
                             "(const date_vector& observation_date, "
                             "const variable_vector2& observation_variable, "
                             "const index_vector3& observation_index, "
                             "const observation_vector3& observation, "
                             "double date, "
                             "variable_vector& aggregated_variable, "
                             "index_vector2& aggregated_index, "
                             "observation_vector2& aggregated_observation)");
    }


    //! Creates a new track.
    /*!
      \return The index of the new track.
    */
    template <class T>
    int ObservationAggregator<T>::CreateTrack()
    {
        Vector<double> track;
        track.PushBack(numeric_limits<double>::min());
        date_.PushBack(track);
        return date_.GetSize() - 1;
    }


    //! Sets the active track to a given track.
    /*!
      \param[in] track the given track.
    */
    template <class T>
    void ObservationAggregator<T>::SetTrack(int track)
    {
        if (track < 0 || track >= date_.GetSize())
            throw WrongIndex("ObservationAggregator<T>::SetTrack(int track)",
                             string("The track should be in [0, ") +
                             to_str(date_.GetSize() - 1) +
                             "], but is equal to " + to_str(track) + ".");
        active_track_index_ = track;
    }


    //! Returns the last date of the current track.
    /*!
      \return The last date of the current track.
    */
    template <class T>
    double ObservationAggregator<T>::LastDate() const
    {
        int Ndate = date_.GetSize(active_track_index_);
        return date_(active_track_index_, Ndate - 1);
    }


    //! Returns the last date of a given track.
    /*!
      \param[in] track a given track.
      \return The last date of the given track.
    */
    template <class T>
    double ObservationAggregator<T>::LastDate(int track) const
    {
        int Ndate = date_.GetSize(track);
        return date_(track, Ndate - 1);
    }


    //! Pushes a date to the current track.
    /*!
      \param[in] date a given date.
    */
    template <class T>
    void ObservationAggregator<T>::PushDate(double date)
    {
        PushDate(date, active_track_index_);
    }



    //! Pushes a date to a given track.
    /*!
      \param[in] track index of a given track.
      \param[in] date a given date.
    */
    template <class T>
    void ObservationAggregator<T>::PushDate(double date, int track)
    {
        int Ndate = date_(track).GetSize();
        date_(track).Resize(Ndate + 1);
        for (int i = Ndate - 1; i > -1; i--)
        {
            date_(track)(i + 1) = date_(track)(i);
            if (date_(track)(i + 1) < date)
            {
                date_(track)(i + 1) = date;
                return;
            }
        }

        date_(0) = date;
    }


    //! Computes the contributions of given observations at a given date.
    /*!
      \param[in] date the given date.
      \param[in] observation_date the dates associated with the given
      observations.
      \param[out] contribution the contributions computed.
    */
    template <class T>
    template <class date_vector>
    void ObservationAggregator<T>
    ::Contribution(double date, const date_vector& observation_date,
                   Vector<double>& contribution) const
    {
        int Ndate = observation_date.GetSize();
        contribution.Reallocate(Ndate);
        for (int i = 0; i< Ndate; i++)
            contribution(i) = Contribution(date - observation_date(i));
    }


    //! Computes the contribution of an observation given a time difference.
    /*!
      \param[in] delta_t the time difference.
    */
    template <class T>
    double ObservationAggregator<T>::Contribution(double delta_t) const
    {
        if (interpolation_type_ == type_step_)
            return 1.;

        if (interpolation_type_ == type_triangle_)
        {
            if (delta_t > 0.)
                return 1. - delta_t / width_left_;
            else
                return 1. - delta_t / width_right_;
        }

        return 0.;
    }


    //! Returns the least interval containing a given value.
    /*! This method returns the closest value \a value_inf among the vector
      \a X elements that is lower than the given value \a value and the
      closest value \a value_sup that is higher.
      \param[in] X a sorted vector.
      \param[in] value a given value.
      \param[out] value_inf lower bound of the interval.
      \param[out] value_sup upper bound of the interval.
    */
    template <class T>
    template <class date_vector>
    void ObservationAggregator<T>::FindInterval(date_vector& X, double value,
                                                double& value_inf,
                                                double& value_sup) const
    {
        int Nx = X.GetSize();
        value_inf = -1.;
        value_sup = -1.;
        for (int i = Nx - 1; i > -1; i--)
            if (X(i) < value)
            {
                value_inf = X(i);
                if (i < Nx - 1)
                    value_sup = X(i + 1);
                return;
            }
    }

}


#define VERDANDI_FILE_OBSERVATION_MANAGER_OSERVATIONAGGREGATOR_CXX
#endif
