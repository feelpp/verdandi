// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet
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



#include <cppunit/extensions/HelperMacros.h>

#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU
#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"
#include "method/BLUE.cxx"
using namespace Verdandi;


class BLUETest: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(BLUETest);
    CPPUNIT_TEST(test_compute_BLUE);
    CPPUNIT_TEST_SUITE_END();

protected:
    int Nx_;
    int Ny_;

public:
    void setUp()
    {
        srand(time(NULL));
    }


    void tearDown()
    {
    }


    // Fills a matrix with ones.
    template <class MatrixType>
    void fill_matrix(int m, int n, MatrixType& M)
    {
        typedef typename MatrixType::value_type T;
        M.Reallocate(m, n);
        M.Fill(T(1));
    }


    // Fills a matrix with ones.
    template <class T, class Prop, class Allocator>
    void fill_matrix(int m, int n, Matrix<T, Prop, RowSparse, Allocator>& M)
    {
        Matrix<T, Prop, ArrayRowSparse, Allocator> Marray(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Marray(i, j) = T(1);
        Copy(Marray, M);
    }


    void test_compute_BLUE()
    {
        int Nx[3] = {10, 10,  1};
        int Ny[3] = { 2, 10,  1};

        for (int i = 0; i < 3; i++)
        {
            Nx_ = Nx[i];
            Ny_ = Ny[i];

            compute_BLUE<Matrix<double, General, RowSparse>,
                Matrix<double, General, RowSparse>,
                Matrix<double, General, RowSparse>,
                Vector<double>,
                Matrix<double, General, RowSparse>,
                Vector<double> >();

            compute_BLUE<Matrix<double>, Matrix<double>, Matrix<double>,
                Vector<double>, Matrix<double>, Vector<double> >();

            compute_BLUE<Matrix<float>, Matrix<float>, Matrix<float>,
                Vector<float>, Matrix<float>, Vector<float> >();
        }
    }


    template <class StateErrorVariance,
              class ObservationOperator, class CrossedMatrix,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void compute_BLUE()
    {
        int i, j;
        {
            StateErrorVariance B(Nx_, Nx_);
            ObservationOperator H(Ny_, Nx_);
            CrossedMatrix tmp;
            ObservationVector y(Ny_);
            ObservationErrorVariance R(Ny_, Ny_);
            StateVector x(Nx_), analysis;

            B.SetIdentity();
            H.SetIdentity();
            R.SetIdentity();
            y.Fill();
            Mlt(3., y);
            x.Fill();

            analysis = x;
            ComputeBLUE_matrix(B, H, tmp, y, R, analysis);

            for (i = 0; i < min(Nx_, Ny_); i++)
                CPPUNIT_ASSERT(analysis(i) == 2. * x(i));
            // No observations for this part of the state.
            for (i = min(Nx_, Ny_); i < Nx_; i++)
                CPPUNIT_ASSERT(analysis(i) == x(i));
        }
        {
            typedef typename StateVector::value_type T;

            StateErrorVariance B(Nx_, Nx_);
            ObservationOperator H;
            CrossedMatrix tmp;
            ObservationVector y(Ny_);
            ObservationErrorVariance R(Ny_, Ny_);
            StateVector x(Nx_), analysis;

            B.SetIdentity();
            fill_matrix(Ny_, Nx_, H);
            R.SetIdentity();
            y.Fill();
            T sum_y = ((Ny_ - 1) * Ny_)/ 2;
            T bg_value = 5;
            x.Fill(bg_value);

            analysis = x;
            ComputeBLUE_matrix(B, H, tmp, y, R, analysis);

            // All elements in the state should have the same value:
            T value = (bg_value + sum_y) / T(1 + Nx_ * Ny_);
            for (i = 0; i < Nx_; i++)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(analysis(i), value,
                                             1.e-6 * value);
        }
    }
};
