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
    QuadraticModel<T>::QuadraticModel(): Delta_t_(1.), date_(0.)
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

        vector<string> state_vector;
        configuration.Set("initial_state", state_vector);
        Nstate_ = int(state_vector.size());
        state_.Reallocate(Nstate_);
        for (int i = 0; i < Nstate_; i++)
            to_num(state_vector[i], state_(i));

        configuration.Set("with_quadratic_term", with_quadratic_term_);
        configuration.Set("with_linear_term", with_linear_term_);
        configuration.Set("with_constant_term", with_constant_term_);

        if (with_quadratic_term_)
        {
            vector<string> Q_vector;
            configuration.Set("quadratic_term", Q_vector);
            if (int(Q_vector.size()) != Nstate_ * Nstate_ * Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"quadratic_term\" has "
                                         + to_str(int(Q_vector.size()))
                                         + " elements, instead of "
                                         + to_str(Nstate_ * Nstate_ * Nstate_)
                                         + ".");
            Q_.resize(Nstate_);
            for (int i = 0; i < Nstate_; i++)
                Q_[i].Reallocate(Nstate_, Nstate_);
            for (int i = 0; i < Nstate_; i++)
                for (int j = 0; j < Nstate_; j++)
                    for (int k = 0; k < Nstate_; k++)
                        to_num(Q_vector[i * Nstate_ * Nstate_
                                        + j * Nstate_ + k], Q_[i](j, k));

            Q_state_.Reallocate(Nstate_);
        }

        if (with_linear_term_)
        {
            vector<string> L_vector;
            configuration.Set("linear_term", L_vector);
            if (int(L_vector.size()) != Nstate_ * Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"linear_term\" has "
                                         + to_str(int(L_vector.size()))
                                         + " elements, instead of "
                                         + to_str(Nstate_ * Nstate_) + ".");
            L_.Reallocate(Nstate_, Nstate_);
            for (int i = 0; i < Nstate_; i++)
                for (int j = 0; j < Nstate_; j++)
                    to_num(L_vector[i * Nstate_ + j], L_(i, j));
        }

        if (with_constant_term_)
        {
            vector<string> b_vector;
            configuration.Set("constant", b_vector);
            if (int(b_vector.size()) != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"constant\" has "
                                         + to_str(int(b_vector.size()))
                                         + " elements, instead of "
                                         + to_str(Nstate_) + ".");
            b_.Reallocate(Nstate_);
            for (int i = 0; i < Nstate_; i++)
                to_num(b_vector[i], b_(i));
        }

        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("initial_date", date_);
        configuration.Set("final_date", final_date_);

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

        date_ += Delta_t_;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool QuadraticModel<T>::HasFinished() const
    {
        return date_ >= final_date_;
    }


    //! Saves the simulated data.
    /*! It saves the state.
     */
    template <class T>
    void QuadraticModel<T>::Save()
    {
        output_saver_.Save(state_, date_, "state");
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


    //! Returns the current date.
    /*!
      \return The current date.
    */
    template <class T>
    double QuadraticModel<T>::GetDate() const
    {
        return date_;
    }


    //! Sets the current date.
    /*!
      \param[in] date the current date.
    */
    template <class T>
    void QuadraticModel<T>::SetDate(double date)
    {
        date_ = date;
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
    ::GetState(typename QuadraticModel<T>::state_vector& state) const
    {
        state = state_;
    }


    //! Sets the controlled state vector.
    /*!
      \param[in] state the new controlled state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::SetState(typename QuadraticModel<T>::state_vector& state)
    {
        state_ = state;
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::GetFullState(typename QuadraticModel<T>::state_vector& state) const
    {
        GetState(state);
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    template <class T>
    void QuadraticModel<T>
    ::SetFullState(const typename QuadraticModel<T>::state_vector& state)
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
