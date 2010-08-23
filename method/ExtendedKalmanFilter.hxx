// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_EXTENDEDKALMANFILTER_HXX


namespace Verdandi
{


    //////////////////////////
    // EXTENDEDKALMANFILTER //
    //////////////////////////


    //! This class implements the extended Kalman filter.
    template <class T, class ClassModel, class ClassObservationManager>
    class ExtendedKalmanFilter: public VerdandiBase
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

    protected:

        /*** Main components ***/

        //! Underlying model.
        ClassModel model_;

        //! Observation manager.
        ClassObservationManager observation_manager_;

        //! Background error covariance matrix (B).
        model_state_error_variance state_error_variance_;

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

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructors and destructor ***/

        ExtendedKalmanFilter(string configuration_file);
        ~ExtendedKalmanFilter();

        /*** Methods ***/

        void Initialize(string configuration_file);

        void InitializeStep();

        void Forward();

        void Analyze();

        void PropagateCovarianceMatrix();
        void PropagateCovarianceMatrix_vector();
        void PropagateCovarianceMatrix_matrix();

        void ComputeBLUE(const observation& innovation, model_state& state);

        bool HasFinished() const;

        // Access methods.
        const ClassModel& GetModel() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_EXTENDEDKALMANFILTER_HXX
#endif
