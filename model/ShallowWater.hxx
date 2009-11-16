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


#ifndef VERDANDI_FILE_SHALLOWWATER_HXX


#include "newran/newran.h"


namespace Verdandi
{


    ////////////////////////
    // SHALLOWWATER MODEL //
    ////////////////////////


    //! This class is a shallow-water model.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class ShallowWater: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef Matrix<T, General, RowSparse> background_error_variance;
        typedef Vector<T> error_covariance_row;
        typedef Vector<T> state_vector;
        typedef Matrix<T, General, RowSparse> crossed_matrix;

    protected:

        //! Water height.
        Matrix<T> h_;
        //! Vertical velocity along x.
        Matrix<T> u_;
        //! Vertical velocity along y.
        Matrix<T> v_;

        //! Water-height flux along x.
        Matrix<T> hf_x_;
        //! Water-height flux along y.
        Matrix<T> hf_y_;
        //! Flux along x of the vertical velocity along x.
        Matrix<T> uf_x_;
        //! Flux along y of the vertical velocity along x.
        Matrix<T> uf_y_;
        //! Flux along x of the vertical velocity along y.
        Matrix<T> vf_x_;
        //! Flux along y of the vertical velocity along y.
        Matrix<T> vf_y_;

        //! First abscissa.
        double x_min_;
        //! First ordinate.
        double y_min_;

        //! Space step along x.
        double Delta_x_;
        //! Space step along y.
        double Delta_y_;

        //! Number of points along x (in the grid for height).
        int Nx_;
        //! Number of points along y (in the grid for height).
        int Ny_;

        //! Time step.
        double Delta_t_;

        //! Number of time steps.
        int Nt_;

        //! Current time step.
        int time_step_;

        //! Gravitational acceleration.
        const double g_;

        /*! \brief Type of boundary condition on the left (0: free, 1: wall,
          2: flow, 3: height). */
        int boundary_condition_left_;
        //! Constant in-flow or height on the left.
        T value_left_;
        //! Amplitude of variations on the left.
        T amplitude_left_;
        //! Frequency of the variations on the left.
        T frequency_left_;

        /*! \brief Type of boundary condition on the right (0: free, 1: wall,
          2: flow, 3: height). */
        int boundary_condition_right_;
        //! Constant in-flow or height on the right.
        T value_right_;
        //! Amplitude of variations on the right.
        T amplitude_right_;
        //! Frequency of the variations on the right.
        T frequency_right_;

        /*! \brief Type of boundary condition on the bottom (0: free, 1: wall,
          2: flow, 3: height). */
        int boundary_condition_bottom_;
        //! Constant in-flow or height on the bottom.
        T value_bottom_;
        //! Amplitude of variations on the bottom.
        T amplitude_bottom_;
        //! Frequency of the variations on the bottom.
        T frequency_bottom_;

        /*! \brief Type of boundary condition on the top (0: free, 1: wall, 2:
          flow, 3: height). */
        int boundary_condition_top_;
        //! Constant in-flow or height on the top.
        T value_top_;
        //! Amplitude of variations on the top.
        T amplitude_top_;
        //! Frequency of the variations on the top.
        T frequency_top_;

        //! Value of the departure in the initial conditions.
        T value_;

        //! Standard deviation of the model error for boundary conditions.
        double model_error_std_bc_;
        //! Standard deviation of the model error for initial conditions.
        double model_error_std_ic_;
        //! Determining the random seed.
        string seed_;
        //! The base uniform random number generator.
        NEWRAN::MotherOfAll* urng_;
        //! Normal random generator.
        NEWRAN::Normal normal_;

        //! Balgovind scale for background covariance.
        double Balgovind_scale_background_;
        //! Background error variance.
        double background_error_variance_value_;
        //! Background error covariance matrix (B).
        background_error_variance background_error_variance_;

        //! Balgovind scale for model covariance.
        double Balgovind_scale_model_;
        //! Model error variance.
        double model_error_variance_;

        //! Number of the row of B currently stored.
        int current_row_;
        //! Number of the column of Q currently stored.
        int current_column_;
        //! Value of the row of B currently stored.
        error_covariance_row error_covariance_row_;

        /*** Experiment settings ***/

        //! Assimilation interval in number of forward integrations.
        int Nt_assimilation_;
        //! Prediction interval in number of forward integrations.
        int Nt_prediction_;
        /*! \brief Flag that indicates whether the positivity of the analyzed
          data is required.
        */
        bool with_positivity_requirement_;

        //! Is there any data not saved yet?
        bool data_to_save_;

        //! Is there any analyzed data not saved yet?
        bool analyzed_data_to_save_;

    public:
        // Constructor and destructor.
        ShallowWater();
        ShallowWater(string configuration_file);
        ~ShallowWater();
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        bool HasFinished() const;
        void StepBack(const state_vector& state);

        // Access methods.
        int GetCurrentDate() const;
        int GetNt() const;
        int GetNx() const;
        int GetNy() const;
        int GetNz() const;
        int GetNs() const;
        int GetXMin() const;
        int GetYMin() const;
        int GetDeltaX() const;
        int GetDeltaY() const;
        int GetNstate() const;
        int GetNtAssimilation() const;
        void GetState(state_vector& state) const;
        void SetState(state_vector& state);
        void GetFullState(state_vector& state) const;
        void SetFullState(const state_vector& state);
        void GetBackgroundErrorCovarianceRow(int row,
                                             error_covariance_row&
                                             error_covariance_row);
        const background_error_variance&
        GetBackgroundErrorVarianceMatrix() const;
        bool IsErrorSparse() const;

        // Access methods useful for saving with the output saver.
        const ShallowWater<T>& GetModel() const;
        bool GetDataToSave() const;
        bool GetAnalyzedDataToSave() const;
        void ClearDataToSave();

    private:
        // Configuration.
        void ReadConfigurationBoundaryCondition(string side,
                                                GetPot& configuration_stream,
                                                int& type, T& value,
                                                T& amplitude, T& frequency);

        //Processing.
        void ComputeGhostCellValue(int type, T value, T amplitude,
                                   T frequency, T h_l, T u_l, T v_l, T& h_r,
                                   T& u_r, T& v_r);
        void ComputeFluxHLL(T h_l, T h_r, T u_l, T u_r, T v_l, T v_r,
                            T& flux_h, T& flux_u, T& flux_v);
        void ComputeFlux(T h, T u, T v, T& flux_h, T& flux_u, T& flux_v);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_SHALLOWWATER_HXX
#endif
