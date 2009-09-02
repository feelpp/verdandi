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


#ifndef VERDANDI_FILE_OUTPUTSAVER_CXX


#include "OutputSaver.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Main constructor.
    /*! Builds the output saver and reads option keys in the configuration
      file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class ClassAssimilationDriver>
    OutputSaver<T, ClassAssimilationDriver>
    ::OutputSaver(string configuration_file, ClassAssimilationDriver& driver)
    {
        GetPot configuration_stream(configuration_file.c_str());


        /***********************
         * Reads configuration *
         ***********************/


        /*** Save options ***/

        configuration_stream.set_prefix("save/");
        output_directory_ = configuration_stream("Output_directory",
                                                 "configuration_error");
        period_save_ = configuration_stream("Period_save", -1);

        /*** Initializations ***/

        h_stream_.open((output_directory_ + "h.bin").c_str());
        h_assimilation_stream_.open((output_directory_
                                     + "h_a.bin").c_str());
        u_stream_.open((output_directory_ + "u.bin").c_str());
        v_stream_.open((output_directory_ + "v.bin").c_str());

        // Saves the initial state.
        Save(driver);

    }


    //! Destructor.
    template <class T, class ClassAssimilationDriver>
    OutputSaver<T, ClassAssimilationDriver>
    ::~OutputSaver()
    {
        h_stream_.close();
        h_assimilation_stream_.close();
        u_stream_.close();
        v_stream_.close();
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the output saver.
    /*!
      \param[in] configuration_file configuration file.
      \param[in] driver the assimilation driver.
    */
    template <class T, class ClassAssimilationDriver>
    void OutputSaver<T, ClassAssimilationDriver>
    ::Initialize(string configuration_file, ClassAssimilationDriver& driver)
    {
        // Saves the state vector with data assimilation.
        Save(driver);
    }


    //! Initializes a step for the output saver.
    template <class T, class ClassAssimilationDriver>
    void OutputSaver<T, ClassAssimilationDriver>
    ::InitializeStep()
    {
    }


    //! Saves the state.
    /*! Saves either the full state or the reduced state after assimilation.
      \param[in] driver the assimilation driver.
    */
    template <class T, class ClassAssimilationDriver>
    void OutputSaver<T, ClassAssimilationDriver>
    ::Save(ClassAssimilationDriver& driver)
    {
        if (driver.GetModel().GetCurrentDate() % period_save_ == 0)
        {

            if (driver.GetDataToSave())
            {

                Vector<T> state_vector;
                driver.GetModel().GetState(state_vector);

                if (driver.GetAnalyzedDataToSave())
                {
                    state_vector.Write(h_assimilation_stream_, false);
                }
                else
                {
                    state_vector.Write(h_stream_, false);
                    state_vector.Write(u_stream_, false);
                    state_vector.Write(v_stream_, false);
                }
                driver.ClearDataToSave();
            }

        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OPTIMALINTERPOLATION_CXX
#endif
