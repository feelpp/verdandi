// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet, Anne Tilloy
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


#ifndef VERDANDI_FILE_METHOD_NEWRANPERTURBATIONMANAGER_HXX

#include "newran.h"

namespace Verdandi
{


    ///////////////////////////////
    // NEWRANPERTURBATIONMANAGER //
    ///////////////////////////////


    //! This class generates random samples using Newran.
    class NewranPerturbationManager: public BasePerturbationManager
    {
    protected:
        //! Uniform random number generator.
        NEWRAN::LGM_mixed* urng_;

        //! Path to the Newran seed directory.
        string seed_path_;

    public:

        /*** Constructors and destructor ***/

        NewranPerturbationManager();
        NewranPerturbationManager(string configuration_file);
        ~NewranPerturbationManager();

        /*** Methods ***/

        void Initialize(string configuration_file);
        void Reinitialize();
        void Finalize();

        using BasePerturbationManager::Normal;
        virtual double Normal(double mean, double variance,
                              Vector<double, VectFull>& parameter);
        virtual void Normal(double mean, double variance,
                            Vector<double, VectFull>& parameter,
                            Vector<double>& sample);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_NEWRANPERTURBATIONMANAGER_HXX
#endif
