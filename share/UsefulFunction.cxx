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


    //! Converts strings to most types.
    /*!
      \param[in] s string to be converted.
      \param[out] out \a s converted to 'T'.
    */
    template <class T>
    void convert(const string& s, T& out)
    {
	istringstream str(s);
	str >> out;
    }


    //! Sets a string.
    /*!
      \param[in] s input string.
      \param[out] out output string, equal to \a s on exit.
    */
    void convert(const string& s, string& out)
    {
	out = s;
    }


    //! Checks whether a string is a number.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is a number, false otherwise.
    */
    bool is_num(const string& str)
    {
	if (str == "")
	    return false;

	bool mant, mant_a, mant_b, exp;
	string::size_type pos;
	string m, e, m_a, m_b;

	pos = str.find_first_of("eE");
	// Mantissa.
	m = str.substr(0, pos);
	// Exponent.
	e = pos == string::npos ? "" : str.substr(pos + 1);

	exp = pos != string::npos;

	pos = m.find_first_of(".");
	// Mantissa in the form: [m_a].[m_b].
	m_a = m.substr(0, pos);
	// Exponent.
	m_b = pos == string::npos ? "" : m.substr(pos + 1);

	mant = m != "" && m != "-" && m != "+";
	mant_a = m_a != "" && m_a != "-" && m_a != "+";
	mant_b = m_b != "";

	return (mant
                && ((mant_a || mant_b)
                    && (!mant_a || is_integer(m_a))
                    && (!mant_b || is_unsigned_integer(m_b)))
                && (!exp || is_integer(e)));
    }


    //! Checks whether a string is an integer.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is an integer, false otherwise.
    */
    bool is_integer(const string& str)
    {
	bool ans;

	ans = (str.size() > 0 && isdigit(str[0]))
	    || (str.size() > 1 && (str[0] == '+' || str[0] == '-'));

	unsigned int i(1);
	while (i < str.size() && ans)
	{
	    ans = ans && isdigit(str[i]);
	    i++;
	}

	return ans;
    }


    //! Checks whether a string is an unsigned integer.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is an unsigned integer, false otherwise.
    */
    bool is_unsigned_integer(const string& str)
    {
	bool ans(str.size() > 0);

	unsigned int i(0);
	while (i < str.size() && ans)
	{
	    ans = ans && isdigit(str[i]);
	    i++;
	}

	return ans;
    }


    //! Trims off a string.
    /*!
      Removes delimiters at each edge of the string.
      \param[in] str string to be trimmed off.
      \param[in] delimiters characters to be removed.
      \return \a str trimmed off.
    */
    string trim(string str, string delimiters)
    {
	string::size_type index_end = str.find_last_not_of(delimiters);
	string::size_type index_beg = str.find_first_not_of(delimiters);

	if (index_beg == string::npos)
	    return "";

	return str.substr(index_beg, index_end - index_beg + 1);
    }


    //! Splits a string.
    /*!
      The string is split according to delimiters and elements are stored
      in the vector 'vect'.
      \param[in] str string to be split.
      \param[out] vect (output) vector containing elements of the string.
      \param[in] delimiters (optional) delimiters. Default: " \n\t".
    */
    template <class T>
    void split(string str, vector<T>& vect, string delimiters)
    {
	vect.clear();

	T tmp;
	string::size_type index_beg, index_end;

	index_beg = str.find_first_not_of(delimiters);

	while (index_beg != string::npos)
	{
	    index_end = str.find_first_of(delimiters, index_beg);
	    convert(str.substr(index_beg, index_end == string::npos ?
			       string::npos : (index_end - index_beg)), tmp);
	    vect.push_back(tmp);
	    index_beg = str.find_first_not_of(delimiters, index_end);
	}
    }


    //! Splits a string.
    /*!
      The string is split according to delimiters.
      \param[in] str string to be split.
      \param[in] delimiters (optional) delimiters. Default: " \n\t".
      \return A vector containing elements of the string.
    */
    vector<string> split(string str, string delimiters)
    {
	vector<string> vect;
	split(str, vect, delimiters);
	return vect;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_USEFULFUNCTION_CXX
#endif
