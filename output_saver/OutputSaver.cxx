// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef VERDANDI_FILE_OUTPUTSAVER_OUTPUT_SAVER_CXX


#include "OutputSaver.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    OutputSaver::OutputSaver():
        save_period_(0), time_tolerance_(0)
    {
    }


    //! Main constructor.
    /*! Reads the configuration.
      \param[in] configuration_file the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    OutputSaver::OutputSaver(string configuration_file, string section_name):
        save_period_(0), time_tolerance_(0)
    {
        Initialize(configuration_file, section_name);
    }


    //! Initializes the output saver with a configuration file.
    /*! Reads the configuration.
      \param[in] configuration_file the configuration file.
      \param[in] section_name the section in \a configuration_file where the
      configuration is to be read.
    */
    void OutputSaver::Initialize(string configuration_file,
                                 string section_name)
    {


        /***********************
         * Reads configuration *
         ***********************/


        GetPot configuration_stream(configuration_file);
        configuration_stream.set_prefix(section_name.c_str());

        configuration_stream.set("Mode", mode_, "", "binary");
        configuration_stream.set("Mode_scalar", mode_scalar_, "", "text");

        string tmp;
        vector<string> variable_vector;
        configuration_stream.set("Variable_list", tmp);
        variable_vector = split(tmp);

        string generic_path;
        configuration_stream.set("File", generic_path);

        vector<string> time_vector;
        configuration_stream.set("Time", tmp, "", "");
        time_vector = split(tmp);
        if ((time_vector.size() >= 1 && time_vector[0] != "step")
            || (time_vector.size() != 0 && time_vector.size() != 2
                && time_vector.size() != 3))
            throw ErrorConfiguration("OutputSaver::OutputSaver(string,"
                                     " string)",
                                     "The variable \"Time\" in file \""
                                     + configuration_file +
                                     "\" cannot be parsed:\n" + tmp);
        if (time_vector.size() >= 2)
            if (!is_num(time_vector[1]))
                throw ErrorConfiguration("OutputSaver::OutputSaver(string,"
                                         " string)",
                                         "The variable \"Time\" in file \""
                                         + configuration_file +
                                         "\" cannot be parsed:\n" + tmp);
            else
                to_num(time_vector[1], save_period_);
        if (time_vector.size() == 3)
            if (!is_num(time_vector[2]))
                throw ErrorConfiguration("OutputSaver::OutputSaver(string,"
                                         " string)",
                                         "The variable \"Time\" in file \""
                                         + configuration_file +
                                         "\" cannot be parsed:\n" + tmp);
            else
                to_num(time_vector[2], time_tolerance_);


        /******************************
         * Initializes variable_list_ *
         ******************************/


        for (unsigned int i = 0; i < variable_vector.size(); i++)
            SetVariable(configuration_stream, generic_path, mode_,
                        variable_vector[i]);
    }


    //! Destructor.
    OutputSaver::~OutputSaver()
    {
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class, that is, "output saver".
    */
    string OutputSaver::GetName() const
    {
        return "output saver";
    }


    /////////////
    // METHODS //
    /////////////


    //! Writes the variable in the file with the chosen format.
    /*! The variable \a variable_name is saved in its dedicated file with its
      saving mode, providing \a time satisfies the proper conditions.
      \param[in] x value of the variable.
      \param[in] time time.
      \param[in] variable_name name of the variable to be saved.
    */
    template <class S>
    void OutputSaver::Save(const S& x, double time, string variable_name)
    {
        if (save_period_ == 0.
            || is_multiple(time, save_period_, time_tolerance_))
            Save(x, variable_name);
    }


    //! Writes the variable in the file with the chosen format.
    /*! The variable \a variable_name is saved in its dedicated file with its
      saving mode.
      \param[in] x value of the variable.
      \param[in] variable_name name of the variable to be saved.
    */
    template <class S>
    void OutputSaver::Save(const S& x, string variable_name)
    {
        map<string, Variable>::iterator im;

        im = variable_list_.find(variable_name);

        if (im == variable_list_.end())
            return;

        Variable& variable = im->second;

        // In case the mode has not been set yet.
        SetVariable<S>(variable);

        if (variable.GetMode() == "text")
            WriteText(x, variable.GetFile());
        else if (variable.GetMode() == "binary")
            WriteBinary(x, variable.GetFile());
    }


    //! Writes \a x in a text file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <class S>
    void OutputSaver::WriteText(const S& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteText(const S& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        file << x << '\n';
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <class S>
    void OutputSaver::WriteBinary(const S& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteBinary(const S& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        x.Write(file);
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <>
    void OutputSaver::WriteBinary(const double& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteBinary(const double& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        file.write(reinterpret_cast<const char*>(&x), sizeof(double));
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <>
    void OutputSaver::WriteBinary(const float& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteBinary(const float& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        file.write(reinterpret_cast<const char*>(&x), sizeof(float));
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <>
    void OutputSaver::WriteBinary(const int& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteBinary(const int& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        file.write(reinterpret_cast<const char*>(&x), sizeof(int));
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
      \warning This method is undefined: it throws an exception in any case.
    */
    template <class T, class Prop, class Allocator>
    void OutputSaver
    ::WriteBinary(const Matrix<T, Prop, RowSparse, Allocator>& x,
                  string file_name) const
    {
        throw ErrorUndefined("WriteBinary(const Matrix<T, Prop, RowSparse, "
                             "Allocator>& x, string file_name)");
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
      \warning This method is undefined: it throws an exception in any case.
    */
    template <class T, class Prop, class Allocator>
    void OutputSaver
    ::WriteBinary(const Matrix<T, Prop, ColSparse, Allocator>& x,
                  string file_name) const
    {
        throw ErrorUndefined("WriteBinary(const Matrix<T, Prop, ColSparse, "
                             "Allocator>& x, string file_name)");
    }


    //! Displays the variables parameters.
    void OutputSaver::DisplayVariableList() const
    {
        map<string, Variable>::const_iterator im;
        for (im = variable_list_.begin(); im != variable_list_.end(); im++)
        {
            cout << "* Variable: " << im->first << "\n";
            im->second.Display();
            cout << endl;
        }
    }


    ////////////////////
    // PRIVATE METHOD //
    ////////////////////


    //! Reads the parameters of the variable in a configuration file.
    /*!
      \param[in] configuration_stream GetPot stream linked to the
      configuration file.
      \param[in] generic_path default output file for all variables.
      \param[in] default_mode default saving format.
      \param[in] variable_name variable name.
    */
    void OutputSaver::SetVariable(GetPot& configuration_stream,
                                  string generic_path,
                                  string default_mode,
                                  string variable_name)
    {
        string current_mode;
        configuration_stream.set("Mode_" + variable_name, current_mode,
                                 "", "");

        string current_path;
        configuration_stream.set("File_" + variable_name, current_path,
                                 "", generic_path);
        current_path = find_replace(current_path, "%{name}", variable_name);

        variable_list_[variable_name] = Variable(current_mode, current_path);
        // If the mode is unknown, the markup "%{extension}" cannot be
        // replaced yet. The markup can only be replaced later, when the type
        // of the variable is known (then the default mode for the given type
        // will be taken). If the mode is already known, the replacement is
        // performed below.
        if (!current_mode.empty())
            SetVariableFile(variable_list_[variable_name]);
    }


    //! Sets the parameters of a variable.
    /*!
      \param[in,out] variable the variable whose parameters are to be
      adjusted.
    */
    template <class S>
    void OutputSaver::SetVariable(Variable& variable)
    {
        if (variable.GetMode().empty())
            // In this case, the mode is not set yet, so the default mode
            // 'mode_' is selected.
        {
            variable.SetMode(mode_);
            SetVariableFile(variable);
        }
    }


    //! Sets the parameters of a variable.
    /*!
      \param[in,out] variable the variable whose parameters are to be
      adjusted.
    */
    template <>
    void OutputSaver::SetVariable<double>(Variable& variable)
    {
        if (variable.GetMode().empty())
            // In this case, the mode is not set yet, so the default mode
            // 'mode_' is selected.
        {
            variable.SetMode(mode_scalar_);
            SetVariableFile(variable);
        }
    }


    //! Sets the parameters of a variable.
    /*!
      \param[in,out] variable the variable whose parameters are to be
      adjusted.
    */
    template <>
    void OutputSaver::SetVariable<float>(Variable& variable)
    {
        if (variable.GetMode().empty())
            // In this case, the mode is not set yet, so the default mode
            // 'mode_' is selected.
        {
            variable.SetMode(mode_scalar_);
            SetVariableFile(variable);
        }
    }


    //! Sets the parameters of a variable.
    /*!
      \param[in,out] variable the variable whose parameters are to be
      adjusted.
    */
    template <>
    void OutputSaver::SetVariable<int>(Variable& variable)
    {
        if (variable.GetMode().empty())
            // In this case, the mode is not set yet, so the default mode
            // 'mode_' is selected.
        {
            variable.SetMode(mode_scalar_);
            SetVariableFile(variable);
        }
    }


    //! Expands the extension in the path to the output file.
    /*! If the file path in \a variable contains the markup "%{extension}",
      this markup is replaced with the extension associated with the saving
      mode.
      \param[in,out] variable the variable whose path may need an expansion of
      the extension.
    */
    void OutputSaver::SetVariableFile(Variable& variable)
    {
        string extension = variable.GetMode();
        if (extension == "binary")
            extension = "bin";
        else if (extension == "text")
            extension = "dat";
        variable.SetFile(find_replace(variable.GetFile(),
                                      "%{extension}", extension));
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OUTPUTSAVER_OUTPUT_SAVER_CXX
#endif
