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


#ifndef VERDANDI_FILE_ERROR_DIAGONALSPARSEMATRIX_HXX


namespace Verdandi
{


    //////////////////////////
    // DIAGONALSPARSEMATRIX //
    //////////////////////////


    //! Diagonal sparse matrix.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class DiagonalSparseMatrix
    {
    protected:

        //! Number of lines (obviously the same as the number of columns).
        int size_;

        //! Value on the diagonal of the matrix.
        T diagonal_value_;

        //! Diagonal sparse matrix structure.
        Matrix<T, General, RowSparse> diagonal_sparse_matrix_;

    public:
        // Constructor and destructor.
        DiagonalSparseMatrix();
        ~DiagonalSparseMatrix();

        // Initialization.
        void Initialize(int size, T diagonal_value);

        // Access.
        int GetSize() const;
        T GetDiagonalValue() const;
        T GetValue(int i, int j) const;
        void GetRow(int row, Vector<T>& matrix_row) const;
        void GetColumn(int column, Vector<T>& matrix_column) const;
        const Matrix<T, General, RowSparse>& GetMatrix() const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_DIAGONALSPARSEMATRIX_HXX
#endif
