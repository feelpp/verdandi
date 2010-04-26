// Copyright (C) 2010 INRIA
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


#ifndef VERDANDI_FILE_METHOD_BLUE_CXX


#include "BLUE.hxx"
#include "seldon/computation/interfaces/direct/SparseSolver.cxx"


namespace Verdandi
{


    //! Computes BLUE using operations on matrices.
    /*! This method is mainly intended for cases where the covariance matrices
      are sparse matrices. Otherwise, the manipulation of the matrices may
      lead to unreasonable memory requirements and to high computational
      costs.
      \param[in] B error variance associated with \a state.
      \param[in] H observation operator.
      \param[in] y the vector of observations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
    */
    template <class StateErrorVariance, class ObservationOperator,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void ComputeBLUE_matrix(const StateErrorVariance& B,
                            const ObservationOperator& H,
                            const ObservationVector& y,
                            const ObservationErrorVariance& R,
                            StateVector& x)
    {
        ComputeBLUE_matrix(B, H, H, y, R, x);
    }


    //! Computes BLUE using operations on matrices.
    /*! This method is mainly intended for cases where the covariance matrices
      are sparse matrices. Otherwise, the manipulation of the matrices may
      lead to unreasonable memory requirements and to high computational
      costs.
      \param[in] B error variance associated with \a state.
      \param[in] H observation operator.
      \param[in] y the vector of observations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
      \param[in] cm this parameter is only used to determine the type of an
      intermediate matrix in the computations. \a cm is not modified nor
      read. Its type CrossedMatrix will be the type of the intermediate matrix
      BH', whose size is Nx times Ny, if Nx is the length of \a x and Ny is
      the length of \a y.
    */
    template <class StateErrorVariance,
              class ObservationOperator, class CrossedMatrix,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void ComputeBLUE_matrix(const StateErrorVariance& B,
                            const ObservationOperator& H,
                            const CrossedMatrix& cm,
                            const ObservationVector& y,
                            const ObservationErrorVariance& R,
                            StateVector& x)
    {
        typedef typename StateVector::value_type T;

        int Ny = y.GetLength();
        int Nx = x.GetLength();

        // Temporary matrix and vector.
        CrossedMatrix working_matrix_xy(Nx, Ny);

        ObservationErrorVariance working_matrix_yy(Ny, Ny);

        // Computes BH'.
        MltAdd(T(1), SeldonNoTrans, B, SeldonTrans, H, T(0),
               working_matrix_xy);

        // Computes HBH'.
        Mlt(H, working_matrix_xy, working_matrix_yy);

        // Computes (HBH' + R).
        Add(T(1), R, working_matrix_yy);

        // Computes inc = (HBH' + R)^{-1} * innovation by solving the linear
        // system (HBH' + R) * inc = innovation.
        ObservationVector innovation = y;
        MltAdd(T(-1), H, x, T(1), innovation);
        Solve(working_matrix_yy, innovation);
        MltAdd(T(1.), working_matrix_xy, innovation, T(1.), x);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BLUE_CXX
#endif
