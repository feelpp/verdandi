// Copyright (C) 2008-2010 INRIA
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


#ifndef VERDANDI_FILE_METHOD_UNSCENTEDKALMANFILTER_HXX


namespace Verdandi
{


    ///////////////////////////
    // UNSCENTEDKALMANFILTER //
    ///////////////////////////


    //! This class implements the unscented Kalman filter.
    template <class T, class ClassModel, class ClassObservationManager>
    class UnscentedKalmanFilter: public VerdandiBase
    {

    public:
        //! Type of a row of the background error variance.
        typedef typename ClassModel::state_error_variance_row
        model_state_error_variance_row;
        //! Type of the model state vector.
        typedef typename ClassModel::state model_state;
        //! Type of the model/observation crossed matrix.
        typedef typename ClassModel::matrix_state_observation
        matrix_state_observation;
        //! Type of the background error variance.
        typedef typename ClassModel::state_error_variance
        model_state_error_variance;
        //! Type of the tangent linear model.
        typedef typename ClassModel::tangent_linear_operator
        model_tangent_linear_operator;
        //! Type of the tangent linear observation operator.
        typedef typename ClassObservationManager
        ::tangent_linear_operator observation_tangent_linear_operator;
        //! Type of a row of the tangent linear observation operator.
        typedef typename ClassObservationManager::tangent_linear_operator_row
        observation_tangent_linear_operator_row;
        //! Type of the observation vector.
        typedef typename ClassObservationManager::observation
        observation;
        //! Type of the sigma point vector.
        typedef Vector<T, VectFull, NewAlloc<T> > sigma_point;
        //! Type of the sigma point collection.
        typedef Vector<sigma_point, Collection> sigma_point_collection;
        //! Type of the state vector collection.
        typedef Vector<model_state, Collection> state_collection;


    protected:

        /*** Main components ***/

        //! Underlying model.
        ClassModel model_;
        //! Observation manager.
        ClassObservationManager observation_manager_;
        //! Background error covariance matrix (B).
        model_state_error_variance background_error_variance_;

        /*** Configuration ***/

        //! Display options.
        map<string, bool> option_display_;
        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Computation mode for BLUE: "vector" or "matrix".
        string blue_computation_;
        //! Computation mode for covariance: "vector" or "matrix".
        string covariance_computation_;

        /*** Sigma-points ***/

        //! Choice of sigma-points.
        string sigma_point_type_;
        //! Collection of sigma-points.
        sigma_point_collection sigma_point_collection_;
        //! Coefficient vector asociated with sigma-points.
        sigma_point alpha_i_;
        //! Boolean to indicate if the coefficients alpha are constants.
        bool alpha_constant_;
        //! alpha.
        T alpha_;

        //! Number of sigma-points.
        int Nsigma_point_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        UnscentedKalmanFilter(string configuration_file);
        ~UnscentedKalmanFilter();

        /*** Methods ***/

        void Initialize(string configuration_file);

        void InitializeStep();

        void Forward();

        void Analyze();

        bool HasFinished() const;

        // Access methods.
        const ClassModel& GetModel() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_UNSCENTEDKALMANFILTER_HXX
#endif