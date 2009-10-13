// Copyright (C) 2008-2009 INRIA
// Author(s): Dominique Chapelle, Philippe Moireau
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


#ifndef VERDANDI_FILE_CLAMPEDBAR_CXX


#include "ClampedBar.hxx"
#include "seldon/SeldonSolver.hxx"

#define DBG_VERBOSE

namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////

    //! Constructor.
    template <class T>
    ClampedBar<T>::ClampedBar()
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    ClampedBar<T>::ClampedBar(string configuration_file)
    {
        //Initialize(configuration_file);
    }


    //! Destructor.
    template <class T>
    ClampedBar<T>::~ClampedBar()
    {
    }

    ////////////////
    // INITIALIZE //
    ////////////////

    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ClampedBar<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        GetPot configuration_stream(configuration_file);

        configuration_stream.set_prefix("domain/");

        configuration_stream.set("bar_length", bar_length_);
        configuration_stream.set("Nx", Nx_);

        configuration_stream.set_prefix("time/");

        configuration_stream.set("Delta_t", Delta_t_);
        configuration_stream.set("time_simu", time_simu_);

        configuration_stream.set_prefix("physics/");

        configuration_stream.set("Young_modulus", Young_modulus_);
        configuration_stream.set("mass_density", mass_density_);

        /*** Allocation ***/

        Delta_x_ = bar_length_/Nx_;
        time_instants_.reserve(floor(time_simu_/Delta_t_));
        time_instants_ = vector<double>(1,0.);
        time_step_ = 0;
        Ndof_ = Nx_+1;

        // Skeleton
        int NvalSkel = 3*(Nx_-1)+4;
        Vector<T> values_0(NvalSkel);
        values_0.Fill(T(0.));
        Vector<int> columns_0(NvalSkel);
        columns_0.Fill();
        Vector<int> rowindex_0(Ndof_+1);
        rowindex_0.Fill();
        columns_0(0) = 0;
        columns_0(1) = 1;
        rowindex_0(0) = 0;
        rowindex_0(1) = 2;
        for (int i=1; i<Nx_; i++)
        {
            columns_0(3*(i-1)+2) = i-1;
            columns_0(3*(i-1)+3) = i;
            columns_0(3*(i-1)+4) = i+1;
            rowindex_0(i+1) = rowindex_0(i)+3;
        }
        columns_0(3*(Nx_-1)+2) = Ndof_-2;
        columns_0(3*(Nx_-1)+3) = Ndof_-1;
        rowindex_0(Ndof_) = 3*(Nx_-1)+4;
        Vector<T> values_1(values_0);
        Vector<int> columns_1(columns_0);
        Vector<int> rowindex_1(rowindex_0);
        Vector<T> values_m(values_0);
        Vector<int> columns_m(columns_0);
        Vector<int> rowindex_m(rowindex_0);
#ifdef DBG_VERBOSE
        DISP(values_0);
        DISP(columns_0);
        DISP(rowindex_0);
#endif
        //Remark: values,rowindex and colums are unlinked after used in
        //SetData therefore, we need one of each for each global matrice
        Newmark_matrix_0_.SetData(Ndof_, Ndof_, values_0,rowindex_0,
                                  columns_0);
        Newmark_matrix_1_.SetData(Ndof_, Ndof_, values_1,rowindex_1,
                                  columns_1);
        Mass_matrix_.SetData(Ndof_, Ndof_, values_m,rowindex_m, columns_m);
#ifdef DBG_VERBOSE
        /* Change values to see the skeleton
	   DISP(Newmark_matrix_0_);
	   DISP(Newmark_matrix_1_);
	   DISP(Mass_matrix_); */
#endif
        // Create local matrix
        Mass_matrix_el_.Reallocate(2,2);
        Stiff_matrix_el_.Reallocate(2,2);
    }

    //! Initializes the first time step for the model.
    template <class T>
    void ClampedBar<T>::InitializeFirstStep()
    {
        /*** Build model ***/
        // Elementary mass matrix construction
        T mass_lin = T(mass_density_*Delta_x_/3);
        Mass_matrix_el_(0,0) = mass_lin;
        Mass_matrix_el_(1,1) = mass_lin;
        Mass_matrix_el_(0,1) = mass_lin/2;
        // Elementary stifness matrix construction
        T stiff_lin = T(Young_modulus_/Delta_x_);
        Stiff_matrix_el_(0,0) = stiff_lin;
        Stiff_matrix_el_(1,1) = stiff_lin;
        Stiff_matrix_el_(0,1) = -stiff_lin;
        // Assembling process
        for (int i=0; i<Nx_; i++)
        {
            //Mass Matrice
            Mass_matrix_.Val(i,i) += Mass_matrix_el_(0,0);
            Mass_matrix_.Val(i+1,i+1) += Mass_matrix_el_(1,1);
            Mass_matrix_.Val(i,i+1) += Mass_matrix_el_(0,1);
            Mass_matrix_.Val(i+1,i) = Mass_matrix_(i,i+1);
            // Newmark's matrice at time n
            Newmark_matrix_0_.Val(i,i) +=
		2.*Mass_matrix_el_(0,0)/(Delta_t_*Delta_t_)+
		-0.5*Stiff_matrix_el_(0,0);
            Newmark_matrix_0_.Val(i+1,i+1) +=
		2.*Mass_matrix_el_(1,1)/(Delta_t_*Delta_t_)+
		-0.5*Stiff_matrix_el_(1,1);
            Newmark_matrix_0_.Val(i,i+1) +=
		2.*Mass_matrix_el_(0,1)/(Delta_t_*Delta_t_)+
		-0.5*Stiff_matrix_el_(0,1);
            Newmark_matrix_0_.Val(i+1,i) = Newmark_matrix_0_(i,i+1);
            //Newmark's matrice at time n+1
            Newmark_matrix_1_.Val(i,i) +=
		2.*Mass_matrix_el_(0,0)/(Delta_t_*Delta_t_)+
		0.5*Stiff_matrix_el_(0,0);
            Newmark_matrix_1_.Val(i+1,i+1) +=
		2.*Mass_matrix_el_(1,1)/(Delta_t_*Delta_t_)+
		0.5*Stiff_matrix_el_(1,1);
            Newmark_matrix_1_.Val(i,i+1) +=
		2.*Mass_matrix_el_(0,1)/(Delta_t_*Delta_t_)+
		0.5*Stiff_matrix_el_(0,1);
            Newmark_matrix_1_.Val(i+1,i) = Newmark_matrix_1_(i,i+1);
        }
#ifdef DBG_VERBOSE
        cout << "Assembled matrices" << endl;
        DISP(Mass_matrix_);
        DISP(Newmark_matrix_1_);
        DISP(Newmark_matrix_0_);
#endif
        // Dirichlet conditions (should be better to read in the skeleton)
        Newmark_matrix_1_.Val(0,0) = 1;
        Newmark_matrix_1_.Val(0,1) = 0;
        Newmark_matrix_1_.Val(1,0) = 0;

        // Vector initialization
        disp_0_.Reallocate(Ndof_);
        disp_1_.Reallocate(Ndof_);
        velo_0_.Reallocate(Ndof_);
        velo_1_.Reallocate(Ndof_);
        force_.Reallocate(Ndof_);
        disp_0_.Fill(T(0.));
        velo_0_.Fill(T(0.));
        disp_1_.Fill(T(0.));
        velo_1_.Fill(T(0.));
        force_.Fill(T(0.));
        // Initial condition
        for (int i=0; i<Ndof_; i++) {
            disp_0_(i) = T(i)/T(Ndof_-1);
        }
    }

    //! Initializes the current time step for the model.
    template <class T>
    void ClampedBar<T>::InitializeStep()
    {
    }

    ////////////////
    // PROCESSING //
    ////////////////

    //! Advances one step forward in time.
    template <class T>
    void ClampedBar<T>::Forward()
    {
        // Update time
        time_instants_.push_back(time_instants_[time_step_]+Delta_t_);
        time_step_ += 1;
        cout << "===================================================" << endl;
        cout << "Iteration " << time_step_ <<
	    " corresponding to time " << time_instants_[time_step_] << endl;
        cout << "--------------------------------------" << endl;

        // Right hand side
        force_.Fill(T(0.));
        MltAdd(2./Delta_t_,Mass_matrix_,velo_0_,1,force_);
        MltAdd(1.,Newmark_matrix_0_,disp_0_,1,force_);
        force_(0) = 0; // Dirichlet conditions

        // initialization of Gmres parameters
        int nb_max_iter = 1000;
        double tolerance = 1e-6;
        Iteration<double> iter(nb_max_iter, tolerance);
        Preconditioner_Base precond; // no preconditioning
        iter.SetRestart(5);
        // Gmres call
        Gmres(Newmark_matrix_1_, disp_1_, force_, precond, iter);

        // velo_1_ = - velo_0_ + 2/Delta_t_(disp_1_-disp_0_);
        velo_1_.Fill(T(0.));
        Add(-1.,velo_0_,velo_1_);
        Add(2./(Delta_t_),disp_1_,velo_1_);
        Add(-2./(Delta_t_),disp_0_,velo_1_);

        // Update
        disp_0_ = disp_1_;
        velo_0_ = velo_1_;
        DISP(disp_0_);

    }

    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool ClampedBar<T>::HasFinished() const
    {
        return time_instants_[time_step_] >= time_simu_;
    }

}

#define VERDANDI_FILE_CLAMPEDBAR_HXX
#endif
