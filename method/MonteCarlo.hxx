// Copyright (C) 2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton, Anne Tilloy
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


#ifndef VERDANDI_FILE_METHOD_MONTECARLO_HXX


#include "BasePerturbationManager.cxx"
#include "NewranPerturbationManager.cxx"


namespace Verdandi
{


    ////////////////
    // MONTECARLO //
    ////////////////


    //! This class performs allows to perform Monte Carlo simulations.
    /*! The class performs a single simulation with perturbed data. In order
      to complete a full Monte Carlo simulation, one has to launch several
      simulations using this class.
    */
    template <class T, class ClassModel>
    class MonteCarlo: public VerdandiBase
    {

    public:
        typedef typename ClassModel::state_vector state_vector;
        typedef typename ClassModel::uncertain_variable uncertain_variable;

    protected:

        //! Underlying model.
        ClassModel model_;

        //! Pertubation managers.
        NewranPerturbationManager perturbation_manager_;

        //! Perturbations vectors.
        vector<uncertain_variable> perturbation_;

        //! Iteration.
        int iteration_;
        //! Date vector.
        Vector<double> date_;

        /*** Configuration ***/

        //! Should the iterations be displayed?
        bool show_iteration_;
        //! Should the current date be displayed?
        bool show_date_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        MonteCarlo(string configuration_file);
        ~MonteCarlo();
        template <class T0, class Storage0, class Allocator0>
        void Clear(Vector<T0, Storage0, Allocator0>& V);
        template <class T0, class Allocator0>
        void Clear(Vector<T0, Collection, Allocator0>& V);

        /*** Main methods ***/

        template <class T0, class Storage0, class Allocator0>
        void SetDimension(Vector<T0, Storage0, Allocator0>& in,
                          Vector<T0, Storage0, Allocator0>& out);

        template <class T0, class Allocator0>
        void SetDimension(Vector<T0, Collection, Allocator0>& in,
                          Vector<T0, Collection, Allocator0>& out);

        template <class T0, class Allocator0>
        void Fill(Vector<T0, Collection, Allocator0>& in, string pdf);

        void Initialize(string configuration_file);
        void InitializeStep();
        void Forward();

        bool HasFinished() const;

        /*** Access methods ***/

        const ClassModel& GetModel() const;
        ClassModel& GetModel();
        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_MONTECARLO_HXX
#endif
