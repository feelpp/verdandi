// Copyright (C) 2009-2010 INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX


namespace Verdandi
{


    /////////////////////
    // QUADRATIC MODEL //
    /////////////////////


    //! This class is a quadratic model.
    /*! The model is defined as \f$\frac{\mathrm{d}x_i}{\mathrm{d}t} = x^T Q_i
      x + L_i x + b_i\f$, where \f$Q_i\f$ is a matrix, \f$L_i\f$ is the
      \f$i\f$-th row of the matrix \f$L\f$ and \f$b\f$ a vector.
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class QuadraticModel: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef Matrix<T> tangent_linear_operator;
        typedef Matrix<T> state_error_variance;
        typedef Vector<T> state_error_variance_row;
        typedef Vector<T> state;
        typedef Matrix<T> matrix_state_observation;
        typedef Matrix<T> error_variance;

    protected:

        //! Dimension of the state.
        int Nstate_;

        //! State vector.
        Vector<T> state_;

        //! Should the quadratic term be applied?
        bool with_quadratic_term_;
        //! Should the linear term be applied?
        bool with_linear_term_;
        //! Should the constant term be applied?
        bool with_constant_term_;

        //! Quadratic terms.
        vector<Matrix<T> > S_;

        //! Matrix that defines the linear part of the model.
        Matrix<T> L_;

        //! Vector that defines the constant part of the model.
        Vector<T> b_;

        //! Time step.
        double Delta_t_;

        //! Final time of the simulation.
        double final_time_;

        //! Current time.
        double time_;

        //! Temporary variable that stores Q times the state vector.
        Vector<T> S_state_;

        /*** Errors ***/

        //! Variance of the model error.
        error_variance Q_;

        //! Variance of the model error in square root form.
        error_variance Q_sqrt_;

        //! Variance of the state error.
        error_variance P_;

        //! Variance of the state error in square root form.
        error_variance P_sqrt_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:
        // Constructors and destructor.
        QuadraticModel();
        QuadraticModel(string configuration_file);
        ~QuadraticModel();
        // Initializations.
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        void ApplyOperator(state& x,
                           bool forward = false, bool preserve_state = true);
        void ApplyTangentLinearOperator(state& x);
        void GetTangentLinearOperator(tangent_linear_operator&) const;
        bool HasFinished() const;
        void Save();

        // Access methods.
        T GetDelta_t() const;
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        void GetState(state& state) const;
        void SetState(state& state);
        void GetFullState(state& state) const;
        void SetFullState(const state& state);

        // Errors.
        error_variance& GetErrorVariance();
        const error_variance& GetErrorVariance() const;
        error_variance& GetErrorVarianceSqrt();
        const error_variance& GetErrorVarianceSqrt() const;
        state_error_variance& GetStateErrorVariance();
        const state_error_variance& GetStateErrorVariance() const;
        void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                      state_error_variance_row);
        state_error_variance& GetStateErrorVarianceSqrt();
        const state_error_variance& GetStateErrorVarianceSqrt() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX
#endif
