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


#ifndef VERDANDI_FILE_SHARE_USEFULFUNCTION_CXX


#include "UsefulFunction.hxx"


namespace Verdandi
{


    //! Returns a value interpolated from a 2D field.
    /*! The output value is produced by bilinear interpolation of the field \a
      input at (\a x, \a y).
      \param x_min abscissa of the first point in \a input.
      \param Delta_x step along x.
      \param y_min ordinate of the first point in \a input.
      \param Delta_y step along y.
      \param input values of the field: \a input(i, j) is the value at
      \f$(x_{min} + i \times \Delta x, y_{min} + j \times \Delta y)\f$.
      \param x abscissa of the point where the field is interpolated.
      \param y ordinate of the point where the field is interpolated.
      \return The value of the field interpolated at (\a x, \a y).
    */
    template <class T, class TM>
    T interpolate(T x_min, T Delta_x, T y_min, T Delta_y,
                  const Matrix<TM>& input, T x, T y)
    {
        int Nx = input.GetN();
        int Ny = input.GetM();

        T distance_x = (x - x_min) / Delta_x;
        int pos_x = int(distance_x);
        int one_pos_x = 1 + pos_x;
        if (pos_x < 0 || pos_x >= Nx - 1)
            throw ErrorConfiguration("interpolate");
        T weight_x = distance_x - T(pos_x);
        T one_weight_x = 1. - weight_x;

        T distance_y = (y - y_min) / Delta_y;
        int pos_y = int(distance_y);
        int one_pos_y = 1 + pos_y;
        if (pos_y < 0 || pos_y >= Ny - 1)
            throw ErrorConfiguration("interpolate");
        T weight_y = distance_y - T(pos_y);
        T one_weight_y = 1. - weight_y;

        return one_weight_y * one_weight_x * input(pos_y, pos_x)
            + one_weight_y * weight_x * input(pos_y, one_pos_x)
            + weight_y * weight_x * input(one_pos_y, one_pos_x)
            + weight_y * one_weight_x * input(one_pos_y, pos_x);
    }


    /*! \brief Returns the position in a multidimensional grid that is
      associated with a global index. */
    /*! A global index gives the position in a grid with a single integer. For
      example, in 2D, the global index of the grid point \f$(i, j)\f$ is \f$i
      \times N + j\f$ if there are N points along the second dimension. This
      function would return \f$(i, j)\f$ in \a position from (\a index =) \f$i
      \times N + j\f$, with \a shape set to \f$(M, N)\f$.
      \param index the global index.
      \param shape dimensions of the grid.
      \param position position in the grid.
    */
    void get_position(int index, const Vector<int>& shape,
                      Vector<int>& position)
    {
        int d;

        int length = shape.GetLength();

        if (index < 0)
            throw ErrorArgument("get_position(int, Vector<int>, Vector<int>)",
                                "Wrong index: " + to_str(index) + ".");
        if (length == 0)
            throw ErrorArgument("get_position(int, Vector<int>, Vector<int>)",
                                "The shape vector is empty.");

        position.Reallocate(length);

        if (length == 1)
            if (index >= shape(0))
                throw ErrorArgument
                    ("get_position(int, Vector<int>, Vector<int>&)",
                     "The shape vector is " + Seldon::to_str(shape)
                     + ", but the index is " + to_str(index) + ".");
            else
            {
                position(0) = index;
                return;
            }

        Vector<int> size(length - 1);
        size(length - 2) = shape(length - 1);
        for (d = length - 3; d >= 0; d--)
            size(d) = size(d + 1) * shape(d + 1);

        for (d = 0; d < length - 1; d++)
        {
            position(d) = index / size(d);
            index = index - position(d) * size(d);
        }
        position(length - 1) = index;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_USEFULFUNCTION_CXX
#endif
