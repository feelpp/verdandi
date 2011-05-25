// Copyright (C) 2010 INRIA
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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_MODEL_MODELTEMPLATE_CXX


#include "ModelTemplate.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    ModelTemplate::ModelTemplate()
    {
    }


    //! Destructor.
    ModelTemplate::~ModelTemplate()
    {
    }


    ////////////////
    // INITIALIZE //
    ////////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    void ModelTemplate::Initialize(string configuration_file)
    {
        throw ErrorUndefined("ModelTemplate"
                             "::Initialize(string configuration_file)");
    }


    //! Initializes the first time step for the model.
    void ModelTemplate::InitializeFirstStep()
    {
        throw ErrorUndefined("ModelTemplate::InitializeFirstStep()");
    }


    //! Initializes the current time step for the model.
    void ModelTemplate::InitializeStep()
    {
        throw ErrorUndefined("ModelTemplate::InitializeStep()");
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    /*! \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f] */
    void ModelTemplate::Forward()
    {
        throw ErrorUndefined("ModelTemplate::Forward()");
    }


    //! Performs one step backward in adjoint model.
    /*!
      \param[in] observation_term \f$ H^T R^{-1}(y - Hx) \f$.
    */
    void ModelTemplate::BackwardAdjoint(state& observation_term)
    {
        throw ErrorUndefined("ModelTemplate::BackwardAdjoint(state&)");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if the simulation is done, false otherwise.
    */
    bool ModelTemplate::HasFinished() const
    {
        throw ErrorUndefined("ModelTemplate::HasFinished() const");
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the model to a given vector.
    /*! The current state of the model is modified.
      \param[in] x a vector.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
    */
    void ModelTemplate::ApplyOperator(state& x,
                                      bool forward, bool preserve_state)
    {
        throw ErrorUndefined("ModelTemplate::ApplyOperator(state& x, "
                             "bool forward, bool preserve_state)");
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in] x a vector.
    */
    void ModelTemplate::ApplyTangentLinearOperator(state& x)
    {
        throw ErrorUndefined("ModelTemplate"
                             "::ApplyTangentLinearOperator(state& x)");
    }


    //! Gets the tangent linear model.
    /*!
      \param[out] A the matrix of the tangent linear model.
    */
    void ModelTemplate
    ::GetTangentLinearOperator(tangent_linear_operator& A) const
    {
        throw ErrorUndefined("ModelTemplate::GetTangentLinearOperator"
                             "(tangent_linear_operator& A) const");
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current time.
    /*!
      \return The current time.
    */
    double ModelTemplate::GetTime() const
    {
        throw ErrorUndefined("ModelTemplate::GetTime() const");
    }


    //! Sets the time of the model to a given time.
    /*!
      \param[in] time a given time.
    */
    void ModelTemplate::SetTime(double time)
    {
        throw ErrorUndefined("ModelTemplate::SetTime(double time)");
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    int ModelTemplate::GetNstate() const
    {
        throw ErrorUndefined("ModelTemplate::GetNstate()");
    }


    //! Returns the size of the full state vector.
    /*!
      \return The size of the full state vector.
    */
    int ModelTemplate::GetNfull_state() const
    {
        throw ErrorUndefined("ModelTemplate::GetNfull_state()");
    }


    //! Provides the state vector.
    /*!
      \param[out] state the reduced state vector.
    */
    void ModelTemplate::GetState(state& state)
    {
        throw ErrorUndefined("ModelTemplate::GetState(state& state)");
    }


    //! Sets the state vector.
    /*! Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
    */
    void ModelTemplate::SetState(state& state)
    {
        throw ErrorUndefined("ModelTemplate"
                             "::SetState(state& state)");
    }


    //! Provides the state lower bound.
    /*!
      \param[out] lower_bound the state lower bound (componentwise).
    */
    void ModelTemplate
    ::GetStateLowerBound(state& lower_bound) const
    {
        throw ErrorUndefined("ModelTemplate"
                             "::GetStateLowerBound(state& lower_bound)");
    }


    //! Provides the state upper bound.
    /*!
      \param[out] upper_bound the state upper bound (componentwise).
    */
    void ModelTemplate
    ::GetStateUpperBound(state& upper_bound) const
    {
        throw ErrorUndefined("ModelTemplate"
                             "::GetStateUpperBound(state& upper_bound)");
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    void ModelTemplate::GetFullState(state& state)
    {
        throw ErrorUndefined("ModelTemplate::GetFullState(state& state)");
    }


    //! Sets the full state vector.
    /*!
      \param[in] state the full state vector.
    */
    void ModelTemplate::SetFullState(state& state)
    {
        throw ErrorUndefined("ModelTemplate::SetFullState(state& state)");
    }


    //! Returns the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    void ModelTemplate::GetAdjointState(state& state_adjoint)
    {
        throw ErrorUndefined("ModelTemplate::GetAdjointState(state& "
                             "state_adjoint)");
    }


    //! Sets the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    void ModelTemplate::SetAdjointState(const state& state_adjoint)
    {
        throw ErrorUndefined("ModelTemplate::SetAdjointState(const state& "
                             "state_adjoint)");
    }


    ////////////
    // ERRORS //
    ////////////


    //! Computes a row of the variance of the state error.
    /*!
      \param[in] row row index.
      \param[out] P_row the row with index \a row in the state error variance.
    */
    void ModelTemplate
    ::GetStateErrorVarianceRow(int row, state_error_variance_row& P_row)
    {
        throw ErrorUndefined("ModelTemplate::GetStateErrorVarianceRow(int"
                             "row, state_error_variance_row&"
                             "state_error_variance_row)");
    }


    //! Returns the state error variance.
    /*!
      \return The state error variance.
    */
    ModelTemplate::state_error_variance&
    ModelTemplate::GetStateErrorVariance()
    {
        throw ErrorUndefined("ModelTemplate::GetStateErrorVariance()");
    }


    /*! Returns a decomposition of the state error covariance matrix (\f$B\f$)
      as a product \f$LUL^T\f$.
    */
    /*!
      \param[out] L the matrix \f$L\f$.
      \param[out] U the matrix \f$U\f$.
    */
    void ModelTemplate::GetStateErrorVarianceSqrt(state_error_variance& L,
                                                  state_error_variance& U)
    {
        throw ErrorUndefined("ModelTemplate::GetStateErrorVarianceSqrt("
                             "state_error_variance& L, "
                             "state_error_variance& U)");
    }


    //! Returns the inverse of the background error variance (\f$B^{-1}\f$).
    /*!
      \return The inverse of the background error variance (\f$B^{-1}\f$).
    */
    const ModelTemplate::state_error_variance&
    ModelTemplate::GetStateErrorVarianceInverse() const
    {
        throw ErrorUndefined("ModelTemplate::GetStateErrorVarianceInverse()");
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    string ModelTemplate::GetName() const
    {
        return "ModelTemplate";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    void ModelTemplate::Message(string message)
    {
        // Put here any processing you need.
    }


}

#define VERDANDI_FILE_MODEL_MODELTEMPLATE_CXX
#endif
