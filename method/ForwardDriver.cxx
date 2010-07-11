// Copyright (C) 2009 INRIA
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


#ifndef VERDANDI_FILE_METHOD_FORWARDDRIVER_CXX


#include "ForwardDriver.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class ClassModel>
    ForwardDriver<ClassModel>::ForwardDriver(string configuration_file):
        model_(configuration_file), iteration_(-1)
    {
        Ops::Ops configuration(configuration_file);

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_,
                                     ClassModel::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ForwardDriver::StaticMessage);

        /*** Display options ***/

        configuration.SetPrefix("forward.display.");
        // Should the iteration be displayed on screen?
        configuration.Set("show_iteration", show_iteration_);
        // Should the date be displayed on screen?
        configuration.Set("show_date", show_date_);

        /*** Ouput saver ***/

        output_saver_.Initialize(configuration_file, "forward.output_saver.");
        output_saver_.Empty("state_forecast");
    }


    //! Destructor.
    template <class ClassModel>
    ForwardDriver<ClassModel>::~ForwardDriver()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the simulation.
    /*! Initializes the model.
      \param[in] configuration_file configuration file to be given to the
      model initialization method.
    */
    template <class ClassModel>
    void ForwardDriver<ClassModel>::Initialize(string configuration_file)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        if (show_date_)
            Logger::StdOut(*this, "Date: " + to_str(model_.GetDate()));
        else
            Logger::Log<-3>(*this,
                            "Date: " + to_str(model_.GetDate()));
        if (show_iteration_)
            Logger::StdOut(*this, "Initialization");
        else
            Logger::Log<-3>(*this, "Initialization");

        model_.Initialize(configuration_file);
        iteration_ = 0;

        MessageHandler::Send(*this, "model", "initial condition");

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes the model before a time step.
    template <class ClassModel>
    void ForwardDriver<ClassModel>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward without optimal interpolation.
    template <class ClassModel>
    void ForwardDriver<ClassModel>::Forward()
    {
        date_.PushBack(model_.GetDate());

        MessageHandler::Send(*this, "all", "::Forward begin");

        if (show_date_)
            Logger::StdOut(*this, "Date: " + to_str(model_.GetDate()));
        else
            Logger::Log<-3>(*this,
                            "Date: " + to_str(model_.GetDate()));
        if (show_iteration_)
            Logger::StdOut(*this, "Iteration " + to_str(iteration_) + " -> "
                           + to_str(iteration_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(iteration_) + " -> "
                            + to_str(iteration_ + 1));

        model_.Forward();
        iteration_++;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class ClassModel>
    bool ForwardDriver<ClassModel>::HasFinished() const
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class ClassModel>
    const ClassModel&
    ForwardDriver<ClassModel>::GetModel() const
    {
        return model_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class ClassModel>
    string ForwardDriver<ClassModel>::GetName() const
    {
        return "ForwardDriver";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class ClassModel>
    void  ForwardDriver<ClassModel>::Message(string message)
    {
        state_vector state;
        if (message.find("forecast") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, model_.GetDate(), "state_forecast");
        }

    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_FORWARDDRIVER_CXX
#endif
