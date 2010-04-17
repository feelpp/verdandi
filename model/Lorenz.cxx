// Copyright (C) 2008-2010 INRIA
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


#ifndef VERDANDI_FILE_MODEL_LORENZ_CXX


#include "Lorenz.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Constructor.
    template <class T>
    Lorenz<T>::Lorenz(): X_(0), Y_(0), Z_(0), Delta_t_(1), date_(0)
    {
    }


    //! Constructor.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    Lorenz<T>::Lorenz(string configuration_file)
    {
        Initialize(configuration_file);
    }


    //! Destructor.
    template <class T>
    Lorenz<T>::~Lorenz()
    {
    }


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    //! Initializes the model.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void Lorenz<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        GetPot configuration_stream(configuration_file, "#", "\n");

        configuration_stream.set_prefix("lorenz/parameter/");

        configuration_stream.set("Prandtl", Pr_);
        configuration_stream.set("Rayleigh", Ra_);
        configuration_stream.set("b", b_);

        configuration_stream.set_prefix("lorenz/initial_condition/");

        configuration_stream.set("X", X_);
        configuration_stream.set("Y", Y_);
        configuration_stream.set("Z", Z_);

        configuration_stream.set_prefix("lorenz/time/");

        configuration_stream.set("Delta_t", Delta_t_);
        configuration_stream.set("Initial_date", date_);
        configuration_stream.set("Final_date", final_date_);

        /*** Output saver ***/

        output_saver_.Initialize(configuration_file, "lorenz/output_saver/");
        output_saver_.Empty("X");
        output_saver_.Empty("Y");
        output_saver_.Empty("Z");
    }


    //! Initializes the current time step for the model.
    template <class T>
    void Lorenz<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    template <class T>
    void Lorenz<T>::Forward()
    {
        X_tmp_ = X_;
        Y_tmp_ = Y_;

        X_ += Delta_t_ * Pr_ * (Y_ - X_);
        Y_ += Delta_t_ * (X_tmp_ * (Ra_ - Z_) - Y_);
        Z_ += Delta_t_ * (X_tmp_ * Y_tmp_ - b_ * Z_);

        date_ += Delta_t_;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool Lorenz<T>::HasFinished() const
    {
        return date_ >= final_date_;
    }


    //! Saves the simulated data.
    /*! It saves the state.
     */
    template <class T>
    void Lorenz<T>::Save()
    {
        output_saver_.Save(X_, date_, "X");
        output_saver_.Save(Y_, date_, "Y");
        output_saver_.Save(Z_, date_, "Z");
    }


    ///////////////////
    // ACCESS METHOD //
    ///////////////////


    //! Returns the value of X.
    /*!
      \return The value of X.
    */
    template <class T>
    T Lorenz<T>::GetX() const
    {
        return X_;
    }


    //! Returns the value of Y.
    /*!
      \return The value of Y.
    */
    template <class T>
    T Lorenz<T>::GetY() const
    {
        return Y_;
    }


    //! Returns the value of Z.
    /*!
      \return The value of Z.
    */
    template <class T>
    T Lorenz<T>::GetZ() const
    {
        return Z_;
    }


    //! Returns the time step.
    /*!
      \return The time step.
    */
    template <class T>
    T Lorenz<T>::GetDelta_t() const
    {
        return Delta_t_;
    }


    //! Returns the current date.
    /*!
      \return The current date.
    */
    template <class T>
    double Lorenz<T>::GetDate() const
    {
        return date_;
    }


    //! Sets the current date.
    /*!
      \param[in] date the current date.
    */
    template <class T>
    void Lorenz<T>::SetDate(double date)
    {
        date_ = date;
    }


    //! Returns the dimension of the state.
    /*!
      \return The dimension of the state, that is, 3.
    */
    template <class T>
    int Lorenz<T>::GetNstate() const
    {
        return 3;
    }


    //! Provides the controlled state vector.
    /*!
      \param[out] state the controlled state vector.
    */
    template <class T>
    void Lorenz<T>::GetState(typename Lorenz<T>::state_vector& state) const
    {
        state.Reallocate(3);
        state(0) = X_;
        state(1) = Y_;
        state(2) = Z_;
    }


    //! Sets the controlled state vector.
    /*!
      \param[in] state the new controlled state vector.
    */
    template <class T>
    void Lorenz<T>::SetState(typename Lorenz<T>::state_vector& state)
    {
        X_ = state(0);
        Y_ = state(1);
        Z_ = state(2);
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    void Lorenz<T>::GetFullState(typename Lorenz<T>::state_vector& state)
        const
    {
        GetState(state);
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    template <class T>
    void Lorenz<T>
    ::SetFullState(const typename Lorenz<T>::state_vector& state)
    {
        SetState(state);
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] error_covariance_row the value of row number \a row.
    */
    template <class T>
    void Lorenz<T>
    ::GetBackgroundErrorCovarianceRow(int row, error_covariance_row&
                                      error_covariance_row)
    {
        throw ErrorUndefined("Lorenz::GetBackgroundErrorCovarianceRow");
    }


    //! Returns the background error covariance matrix (B) if available.
    /*! Returns the background error covariance matrix (B) if available,
      raises an exception otherwise.
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename Lorenz<T>::background_error_variance& Lorenz<T>
    ::GetBackgroundErrorVarianceMatrix() const
    {
        throw ErrorUndefined("Lorenz::GetBackgroundErrorVarianceMatrix");
    }


    //! Checks if the error covariance matrix is sparse.
    /*!
      \return True if there is a sparse error matrix, false otherwise.
    */
    template <class T>
    bool Lorenz<T>::IsErrorSparse() const
    {
        throw ErrorUndefined("Lorenz::IsErrorSparse");
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string Lorenz<T>::GetName() const
    {
        return "Lorenz";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void Lorenz<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_LORENZ_CXX
#endif
