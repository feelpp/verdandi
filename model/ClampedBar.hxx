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


#ifndef VERDANDI_FILE_CLAMPEDBAR_HXX


#include "seldon/SeldonSolver.hxx"

#include "OutputSaver.cxx"




namespace Verdandi
{

    //////////////////////
    // CLAMPEDBAR MODEL //
    //////////////////////


    //! This class is a clamped-bar model.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class ClampedBar: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

    protected:
        //! Bar length.
        double bar_length_;
        //! Space step along x.
        double Delta_x_;
        //! Number of elements along x.
        int Nx_;
        //! Number of degrees of freedom (dofs)
        int Ndof_;

        //! Simulation duration.
        double time_simu_;
        //! Time step.
        double Delta_t_;
        //! Current time step.
        int time_step_;
        //! Time steps
        vector<double> time_instants_;

        //! Mass parameter.
        double mass_density_;
        //! Young's Modulus.
        double Young_modulus_;


        //! FEM Vector
        Vector<T> disp_0_;
        Vector<T> velo_0_;
        Vector<T> disp_1_;
        Vector<T> velo_1_;
        Vector<T> force_;

        //! Mass FEM matrix.
        Matrix<T, General, RowMajor> Mass_matrix_el_;
        //! Stiffness FEM matrix.
        Matrix<T, General, RowMajor> Stiff_matrix_el_;

        //! Newmark Global FEM matrix
        Matrix<T, Symmetric, RowSymSparse> Mass_matrix_;
        Matrix<T, Symmetric, RowSymSparse> Newmark_matrix_0_;
        Matrix<T, Symmetric, RowSymSparse> Newmark_matrix_1_;



#if defined(SELDON_WITH_UMFPACK) && defined(VERDANDI_WITH_DIRECT_SOLVER)
        MatrixUmfPack<T> mat_lu;
#endif
#if defined(SELDON_WITH_SUPERLU) && defined(VERDANDI_WITH_DIRECT_SOLVER)
        MatrixSuperLU<T> mat_lu;
#endif
#if defined(SELDON_WITH_MUMPS) && defined(VERDANDI_WITH_DIRECT_SOLVER)
        MatrixMumps<T> mat_lu;
#endif




    public:
        // Constructor and destructor.
        ClampedBar();
        ClampedBar(string configuration_file);
        ~ClampedBar();
        void Initialize(string configuration_file);
        void InitializeFirstStep();
        void InitializeStep();

        // Processing.
        void Forward();
        bool HasFinished() const;

        string GetName() const;


    };


} // namespace Verdandi.


#define VERDANDI_FILE_CLAMPEDBAR_HXX
#endif
