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


#ifndef VERDANDI_FILE_OPTIMALINTERPOLATION_HXX


#include <map>

#include "ShallowWater.hxx"
#include "GridToNetworkObservationManager.hxx"
#include "Seldon.hxx"
#include "Talos.hxx"


namespace Verdandi
{


    using namespace std;
    using namespace Seldon;
    using namespace Talos;


    //////////////////////////
    // OPTIMALINTERPOLATION //
    //////////////////////////


    template <class T, class ClassModel, class ClassObservationManager>
    class OptimalInterpolation
    {

    protected:

        /*** Main components ***/

        //! Underlying model.
        ClassModel model_;

        //! Observation manager.
        ClassObservationManager observation_manager_;

        /*** Configuration ***/

        //! Display options.
        map<string, bool> option_display_;

        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;

        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;

        //! Is there any data not saved yet?
        bool data_to_save_;

        //! Is there any analyzed data not saved yet?
        bool analyzed_data_to_save_;

    public:

        /*** Constructors and destructor ***/

        OptimalInterpolation();
        OptimalInterpolation(string configuration_file);
        ~OptimalInterpolation();

        /*** Methods ***/

        void Initialize(string configuration_file);

        void InitializeStep();

        void Forward();

        void Analyze();

        void ComputeBLUE(Vector<T>& state_vector);

        bool HasFinished() const;

        // Access methods.
        const ClassModel& GetModel() const;
        bool GetDataToSave() const;
        bool GetAnalyzedDataToSave() const;
        void ClearDataToSave();

    };


} // namespace Verdandi.


#define VERDANDI_FILE_OPTIMALINTERPOLATION_HXX
#endif
