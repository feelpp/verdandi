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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_OBSERVATIONAGGREGATOR_HXX


namespace Verdandi
{


    ///////////////////////////
    // OSERVATION AGGREGATOR //
    ///////////////////////////


    //! Observation manager which computes an aggregated observation vector.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class ObservationAggregator
    {
    private:

        /*** Tracks ***/

        //! Vector of dates.
        Vector2<double> date_;
        //! Index of the active track.
        int active_track_index_;

        /*** Aggregation ***/

        enum interpolation_type {type_step_, type_triangle_,
                                 type_interpolation_};
        //! Interpolation function.
        interpolation_type interpolation_type_;

        enum width_property {width_constant_, width_per_observation_};
        //! Width property.
        width_property width_property_;

        //! Step width or constant triangle left width.
        double width_left_;
        //! Step width or constant triangle rigth width.
        double width_right_;

        //! If the triangle widths are not constant, one should define an
        //! observation interval. We assume that the observations outside
        //! this interval have no contribution.
        double width_left_upper_bound_;
        double width_right_upper_bound_;

        //! The maximal contribution of each observation.
        bool discard_observation_;


    public:

        /*** Constructors and destructor ***/

        ObservationAggregator();
        ObservationAggregator(string configuration_file);
        ~ObservationAggregator();

        void Initialize(string configuration_file);

        void GetContributionInterval(double date, double& date_inf,
                                     double& date_sup, int& selection_policy)
            const;

        template <class date_vector, class observation_vector2,
                  class observation_vector>
        void Aggregate(const date_vector& observation_date,
                       const Vector<double>& contribution,
                       const observation_vector2& observation,
                       double date,
                       observation_vector& aggregated_observation);
        template <class date_vector, class variable_vector2,
                  class observation_vector3,
                  class variable_vector, class observation_vector2>
        void Aggregate(const date_vector& observation_date,
                       const Vector<double>& contribution,
                       const variable_vector2& observation_variable,
                       const observation_vector3& observation,
                       double date,
                       variable_vector& aggregated_variable,
                       observation_vector2& aggregated_observation);
        template <class date_vector, class variable_vector2,
                  class index_vector3, class observation_vector3,
                  class variable_vector, class index_vector2,
                  class observation_vector2>
        void Aggregate(const date_vector& observation_date,
                       const Vector<double>& contribution,
                       const variable_vector2& observation_variable,
                       const index_vector3& observation_index,
                       const observation_vector3& observation,
                       double date,
                       variable_vector& aggregated_variable,
                       index_vector2& aggregated_index,
                       observation_vector2& aggregated_observation);

        /*** Tracks management ***/

        int CreateTrack();
        void SetTrack(int track);
        double LastDate() const;
        double LastDate(int track) const;
        void PushDate(double date);
        void PushDate(double date, int track);

        /*** Contributions management ***/

        template <class date_vector>
        void Contribution(double date, const date_vector& observation_date,
                          Vector<double>& contribution);
        template <class date_vector>
        void Contribution(double date, const date_vector& observation_date,
                          Vector<double>& width_left,
                          Vector<double>& width_right,
                          Vector<double>& contribution);
        double Contribution(double delta_t) const;

        template <class date_vector>
        void GetValueIndex(date_vector& X, double value, int& index_inf,
                           int& index_sup) const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_OBSERVATIONAGGREGATOR_HXX
#endif
