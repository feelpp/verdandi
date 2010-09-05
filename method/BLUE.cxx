// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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
      \param[in] y the vector of observations or innovations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
      \param[in] is_y_innovation Boolean to indicate if the parameter \a y is
      a vector of observations or innovations.
      \param[in] compute_variance Boolean to indicate if the covariance matrix
      has to be updated.
    */
    template <class StateErrorVariance, class ObservationOperator,
              class Observation, class ObservationErrorVariance,
              class State>
    void ComputeBLUE_matrix(StateErrorVariance& B,
                            const ObservationOperator& H,
                            const Observation& y,
                            const ObservationErrorVariance& R,
                            State& x,
                            bool is_y_innovation,
                            bool compute_variance)
    {
        ComputeBLUE_matrix(B, H, H, y, R, x, is_y_innovation,
                           compute_variance);
    }


    //! Computes BLUE using operations on matrices.
    /*! This method is mainly intended for cases where the covariance matrices
      are sparse matrices. Otherwise, the manipulation of the matrices may
      lead to unreasonable memory requirements and to high computational
      costs.
      \param[in] B error variance associated with \a state.
      \param[in] H observation operator.
      \param[in] cm this parameter is only used to determine the type of an
      intermediate matrix in the computations. \a cm is not modified nor
      read. Its type CrossedMatrix will be the type of the intermediate matrix
      BH', whose size is Nx times Ny, if Nx is the length of \a x and Ny is
      the length of \a y.
      \param[in] y the vector of observations or innovations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
      \param[in] is_y_innovation Boolean to indicate if the parameter \a y is
      a vector of observations or innovations.
      \param[in] compute_variance Boolean to indicate if the covariance matrix
      has to be updated.
    */
    template <class StateErrorVariance,
              class ObservationOperator, class MatrixStateObservation,
              class Observation, class ObservationErrorVariance,
              class State>
    void ComputeBLUE_matrix(StateErrorVariance& B,
                            const ObservationOperator& H,
                            const MatrixStateObservation& cm,
                            const Observation& y,
                            const ObservationErrorVariance& R,
                            State& x,
                            bool is_y_innovation,
                            bool compute_variance)
    {
        typedef typename State::value_type T;

        int Ny = y.GetLength();
        int Nx = x.GetLength();

        // Temporary matrices.
        MatrixStateObservation working_matrix_xy(Nx, Ny);

        ObservationErrorVariance working_matrix_yy(Ny, Ny);

        // Computes BH'.
        MltAdd(T(1), SeldonNoTrans, B, SeldonTrans, H, T(0),
               working_matrix_xy);

        // Computes HBH'.
        Mlt(H, working_matrix_xy, working_matrix_yy);

        // Computes (HBH' + R).
        Add(T(1), R, working_matrix_yy);

        // Innovation.
        Observation innovation;
        if (is_y_innovation)
            innovation.SetData(y);
        else
        {
            innovation = y;
            MltAdd(T(-1), H, x, T(1), innovation);
        }

        if (!compute_variance)
        {
            // Computes inc = (HBH' + R)^{-1} * innovation by solving the
            // linear system (HBH' + R) * inc = innovation.
            GetAndSolveLU(working_matrix_yy, innovation);
            MltAdd(T(1), working_matrix_xy, innovation, T(1), x);
        }
        else
        {
            // Kalman Gain.
            MatrixStateObservation K(Nx, Ny);
            K.Fill(T(0));

            GetInverse(working_matrix_yy);

            MltAdd(T(1), working_matrix_xy, working_matrix_yy, T(0), K);

            MltAdd(T(1), K, innovation, T(1), x);

            MltAdd(T(-1), SeldonNoTrans, K, SeldonTrans,
                   working_matrix_xy, T(1), B);
        }

        if (is_y_innovation)
            innovation.Nullify();
    }

} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BLUE_CXX
#endif
