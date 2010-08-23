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


#ifndef VERDANDI_FILE_MODEL_QUADRATICMODEL_CXX


#include "QuadraticModel.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Constructor.
    template <class T>
    QuadraticModel<T>::QuadraticModel(): Delta_t_(1.), time_(0.)
    {
    }


    //! Constructor.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    QuadraticModel<T>::QuadraticModel(string configuration_file)
    {
        Initialize(configuration_file);
    }


    //! Destructor.
    template <class T>
    QuadraticModel<T>::~QuadraticModel()
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
    void QuadraticModel<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        Ops configuration(configuration_file);

        configuration.SetPrefix("quadratic_model.definition.");

        configuration.Set("initial_state", state_);
        Nstate_ = state_.GetLength();

        configuration.Set("with_quadratic_term", with_quadratic_term_);
        configuration.Set("with_linear_term", with_linear_term_);
        configuration.Set("with_constant_term", with_constant_term_);

        if (with_quadratic_term_)
        {
            Q_state_.Reallocate(Nstate_);
            Q_.resize(Nstate_);
            for (int i = 0; i < Nstate_; i++)
                Q_[i].Reallocate(Nstate_, Nstate_);
            configuration.Set("quadratic_term", Q_);
            if (int(Q_.size()) != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"quadratic_term\" has "
                                         + to_str(int(Q_.size()))
                                         + " elements.");
            for (int i = 0; i < Nstate_; i++)
                if (Q_[i].GetM() != Nstate_ || Q_[i].GetN() != Nstate_)
                    throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                             "The initial state has "
                                             + to_str(Nstate_) + " elements, "
                                             "but the matrix " + to_str(i)
                                             + " of \"quadratic_term\" has "
                                             + to_str(int(Q_[i].GetM()))
                                             + " rows and "
                                             + to_str(int(Q_[i].GetN()))
                                             + " columns.");
        }

        if (with_linear_term_)
        {
            L_.Reallocate(Nstate_, Nstate_);
            configuration.Set("linear_term", L_);
            if (L_.GetM() != Nstate_ || L_.GetN() != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"linear_term\" is a "
                                         + to_str(int(L_.GetM())) + " x "
                                         + to_str(int(L_.GetN()))
                                         + " matrix.");
        }

        if (with_constant_term_)
        {
            configuration.Set("constant", b_);
            if (b_.GetLength() != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"constant\" has "
                                         + to_str(b_.GetLength())
                                         + " elements.");
        }

        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("initial_time", time_);
        configuration.Set("final_time", final_time_);

        /*** Output saver ***/

        output_saver_.Initialize(configuration_file,
                                 "quadratic_model.output_saver.");
        if (with_quadratic_term_)
        {
            output_saver_.Empty("Q");
            for (int i = 0; i < Nstate_; i++)
                output_saver_.Save(Q_[i], "Q");
        }
        if (with_linear_term_)
        {
            output_saver_.Empty("L");
            output_saver_.Save(L_, "L");
        }
        if (with_constant_term_)
        {
            output_saver_.Empty("b");
            output_saver_.Save(b_, "b");
        }
        output_saver_.Empty("state");
    }


    //! Initializes the current time step for the model.
    template <class T>
    void QuadraticModel<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    template <class T>
    void QuadraticModel<T>::Forward()
    {
        if (with_quadratic_term_)
            for (int i = 0; i < Nstate_; i++)
            {
                MltAdd(Delta_t_, Q_[i], state_, T(0), Q_state_);
                state_(i) += DotProd(Q_state_, state_);
            }
        if (with_linear_term_)
            MltAdd(Delta_t_, L_, state_, T(1), state_);
        if (with_constant_term_)
            Add(Delta_t_, b_, state_);

        time_ += Delta_t_;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool QuadraticModel<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Saves the simulated data.
    /*! It saves the state.
     */
    template <class T>
    void QuadraticModel<T>::Save()
    {
        output_saver_.Save(state_, time_, "state");
    }


    ///////////////////
    // ACCESS METHOD //
    ///////////////////


    //! Returns the time step.
    /*!
      \return The time step.
    */
    template <class T>
    T QuadraticModel<T>::GetDelta_t() const
    {
        return Delta_t_;
    }


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T>
    double QuadraticModel<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the current time.
    /*!
      \param[in] time the current time.
    */
    template <class T>
    void QuadraticModel<T>::SetTime(double time)
    {
        time_ = time;
    }


    //! Returns the dimension of the state.
    /*!
      \return The dimension of the state.
    */
    template <class T>
    int QuadraticModel<T>::GetNstate() const
    {
        return Nstate_;
    }


    //! Provides the controlled state vector.
    /*!
      \param[out] state the controlled state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::GetState(typename QuadraticModel<T>::state& state) const
    {
        state = state_;
    }


    //! Sets the controlled state vector.
    /*!
      \param[in] state the new controlled state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::SetState(typename QuadraticModel<T>::state& state)
    {
        state_ = state;
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::GetFullState(typename QuadraticModel<T>::state& state) const
    {
        GetState(state);
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::SetFullState(const typename QuadraticModel<T>::state& state)
    {
        SetState(state);
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string QuadraticModel<T>::GetName() const
    {
        return "QuadraticModel";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void QuadraticModel<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_QUADRATICMODEL_CXX
#endif
