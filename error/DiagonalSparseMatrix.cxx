// Copyright (C) 2009 INRIA
// Author(s): Claire Mouton, Vivien Mallet
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


#ifndef VERDANDI_FILE_ERROR_DIAGONALSPARSEMATRIX_CXX


#include "DiagonalSparseMatrix.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    template <class T>
    DiagonalSparseMatrix<T>
    ::DiagonalSparseMatrix()
    {
    }


    //! Destructor.
    template <class T>
    DiagonalSparseMatrix<T>
    ::~DiagonalSparseMatrix()
    {
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the diagonal sparse matrix.
    /*!
      param[in] size size of the matrix.
      param[in] diagonal_value value on the diagonal of the matrix.
    */
    template <class T>
    void DiagonalSparseMatrix<T>
    ::Initialize(int size, T diagonal_value)
    {
        size_ = size;
        diagonal_value_ = diagonal_value;

        /*** Building the diagonal sparse matrix ***/

        Vector<int> column(size_);
        Vector<int> pointer(size_ + 1);
        Vector<T> value(size_);

        value.Fill(diagonal_value_);

        for (int i = 0; i < size_; i++)
        {
            column(i) = i;
            pointer(i) = i;
        }
        pointer(size_) = size_;

        diagonal_sparse_matrix_.SetData(size_, size_, value, pointer, column);
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the number of lines.
    /*!
      \return Number of lines.
    */
    template <class T>
    int DiagonalSparseMatrix<T>::GetSize() const
    {
        return size_;
    }


    //! Returns the value on the diagonal of the matrix.
    /*!
      \return Value on the diagonal of the matrix.
    */
    template <class T>
    T DiagonalSparseMatrix<T>::GetDiagonalValue() const
    {
        return diagonal_value_;
    }


    //! Value of an element of the diagonal sparse matrix.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the diagonal sparse matrix.
    */
    template <class T>
    T DiagonalSparseMatrix<T>::GetValue(int i, int j) const
    {
        if (i == j)
            return diagonal_value_;
        else
            return T(0);
    }


    //! One row of the diagonal sparse matrix structure.
    /*!
      \param[in] row row index.
      \param[out] matrix_row the row \a row of the diagonal sparse matrix.
    */
    template <class T>
    void DiagonalSparseMatrix<T>::GetRow(int row, Vector<T>& matrix_row) const
    {
        matrix_row.Reallocate(size_);
        matrix_row.Zero();
        matrix_row(row) = diagonal_value_;
    }


    //! One column of the diagonal sparse matrix structure.
    /*!
      \param[in] column column index.
      \param[out] matrix_column the column \a column of the diagonal sparse
      matrix.
    */
    template <class T>
    void DiagonalSparseMatrix<T>
    ::GetColumn(int column, Vector<T>& matrix_column) const
    {
        GetRow(column, matrix_column);
    }


    //! Diagonal sparse matrix structure.
    /*!
      \return The diagonal sparse matrix structure.
    */
    template <class T>
    const Matrix<T, General, RowSparse>& DiagonalSparseMatrix<T>
    ::GetMatrix() const
    {
        return diagonal_sparse_matrix_;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_DIAGONALSPARSEMATRIX_CXX
#endif
