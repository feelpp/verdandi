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


#ifndef VERDANDI_FILE_MODEL_SHALLOWWATER_HXX


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
        //! The numerical type (e.g., double).
        typedef T value_type;
        //! Pointer to the numerical type.
        typedef T* pointer;
        //! Const pointer to the numerical type.
        typedef const T* const_pointer;
        //! Reference to the numerical type.
        typedef T& reference;
        //! Const reference to the numerical type.
        typedef const T& const_reference;
        //! Type of the background error covariance matrix.
        typedef Matrix<T, General, RowSparse> state_error_variance;
        //! Type of a row of the background error variance.
        typedef Vector<T> state_error_variance_row;
        //! Type of the model state vector.
        typedef Vector<T> state;
        //! Type of the model/observation crossed matrix.
        typedef Matrix<T, General, RowSparse> matrix_state_observation;

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
        //! Current time.
        double time_;
        //! Simulation duration.
        double final_time_;

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
        double state_error_variance_value_;
        //! Background error covariance matrix (B).
        state_error_variance state_error_variance_;

        //! Balgovind scale for model covariance.
        double Balgovind_scale_model_;
        //! Model error variance.
        double model_error_variance_;

        //! Number of the row of B currently stored.
        int current_row_;
        //! Number of the column of Q currently stored.
        int current_column_;
        //! Value of the row of B currently stored.
        state_error_variance_row state_error_variance_row_;

        /*** Experiment settings ***/

        /*! \brief Flag that indicates whether the positivity of the analyzed
          data is required.
        */
        bool with_positivity_requirement_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

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
        void StepBack(const state& state);
        void Save();

        // Access methods.
        double GetTime() const;
        int GetNx() const;
        int GetNy() const;
        int GetNz() const;
        int GetXMin() const;
        int GetYMin() const;
        int GetDeltaX() const;
        int GetDeltaY() const;
        int GetNstate() const;
        void GetState(state& state) const;
        void SetState(state& state);
        void GetFullState(state& state) const;
        void SetFullState(const state& state);
        void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                      state_error_variance_row);
        const state_error_variance& GetStateErrorVariance() const;
        bool IsErrorSparse() const;

        const ShallowWater<T>& GetModel() const;

        string GetName() const;
        void Message(string message);

    private:
        // Configuration.
        void
        ReadConfigurationBoundaryCondition(string side,
                                           Ops& configuration,
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


#define VERDANDI_FILE_MODEL_SHALLOWWATER_HXX
#endif
