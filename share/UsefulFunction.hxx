// Copyright (C) 2008, INRIA
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


#ifndef VERDANDI_FILE_SHARE_USEFULFUNCTION_HXX

#include <sstream>
#include <string>
#include <iostream>

namespace Verdandi
{


    template <class T, class TM>
    T interpolate(T x_min, T Delta_x, T y_min, T Delta_y,
                  const Matrix<TM>& input, T x, T y);


    void get_position(int index, const Vector<int>& shape,
                      Vector<int>& position);

    /*** From Talos library ***/

    bool is_num(const string& s);
    bool is_integer(const string& s);
    bool is_unsigned_integer(const string& s);
    string trim(string str, string delimiters = " \n\t");
    template <class T>
    void split(string str, vector<T>& vect, string delimiters = " \n\t");
    vector<string> split(string str, string delimiters = " \n\t");
    string find_replace(string str, string old_str, string new_str);
    template<class T>
    bool is_equal(T x, T y, T epsilon = 1.e-6);
    template<class T>
    bool is_multiple(T x, T d, T epsilon = 1.e-6);
    string upper_case(string str);

    /*** Builds a diagonal sparse matrix ***/

    template <class T>
    void build_diagonal_sparse_matrix(int size, T diagonal_value,
                                      Matrix<T, General, RowSparse>& matrix);


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_USEFULFUNCTION_HXX
#endif
