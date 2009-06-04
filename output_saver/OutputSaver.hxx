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


#ifndef VERDANDI_FILE_OUTPUTSAVER_HXX


#include "Seldon.hxx"
#include "Talos.hxx"


namespace Verdandi
{


    using namespace std;
    using namespace Seldon;
    using namespace Talos;


    /////////////////
    // OUTPUTSAVER //
    /////////////////


    template <class T, class ClassAssimilationDriver>
    class OutputSaver
    {

    protected:

        /*** Main components ***/

        //! Output directory.
        string output_directory_;
        //! Number of time steps between two saves.
        int period_save_;
        //! Output ofstream for h.
        ofstream h_stream_;
        //! Output ofstream for h after assimilation.
        ofstream h_assimilation_stream_;
        //! Output ofstream for u.
        ofstream u_stream_;
        //! Output ofstream for v.
        ofstream v_stream_;

    public:

        /*** Constructors and destructor ***/

        OutputSaver();
        OutputSaver(string configuration_file,
                    ClassAssimilationDriver& driver);
        ~OutputSaver();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        ClassAssimilationDriver& driver);

        void InitializeStep();

        void Save(ClassAssimilationDriver& driver);

    };


} // namespace Verdandi.


#define VERDANDI_FILE_OUTPUTSAVER_HXX
#endif
