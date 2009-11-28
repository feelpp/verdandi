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


#ifndef VERDANDI_FILE_OUTPUTSAVER_OUTPUTSAVER_HXX


#include "Variable.cxx"

#include <iostream>
#include <fstream>
#include <string>


namespace Verdandi
{


    /////////////////
    // OUTPUTSAVER //
    /////////////////



    class OutputSaver
    {

    private:

        //! Time interval between two saves.
        double save_period_;
        //! Tolerance on the time interval.
        double time_tolerance_;

        //! Default mode.
        string mode_;

        //! Default mode for scalar variables.
        string mode_scalar_;

        //! Stores the variables properties.
        map<string, Variable> variable_list_;

    public:

        /*** Constructors and destructor ***/

        OutputSaver();
        OutputSaver(string configuration_file, string method_name);
        void Initialize(string configuration_file, string method_name);

        ~OutputSaver();

        string GetName() const;

        /*** Methods ***/

        template <class S>
        void Save(const S& x, double time, string variable_name);

        template <class S>
        void Save(const S& x, string variable_name);

        template <class S>
        void WriteText(const S& x, string file_name) const;

        template <class S>
        void WriteBinary(const S& x, string file_name) const;

        template <class T, class Prop, class Allocator>
        void WriteBinary(const Matrix<T, Prop, RowSparse, Allocator>& x,
                         string file_name) const;

        template <class T, class Prop, class Allocator>
        void WriteBinary(const Matrix<T, Prop, ColSparse, Allocator>& x,
                         string file_name) const;

        void DisplayVariableList() const;

    private:

        void SetVariable(GetPot& configuration_stream,
                         string generic_path,
                         string default_mode,
                         string variable_name);
        template <class S>
        void SetVariable(Variable& variable);
        void SetVariableFile(Variable& variable);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_OUTPUTSAVER_OUTPUTSAVER_HXX
#endif

