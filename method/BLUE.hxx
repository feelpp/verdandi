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


#ifndef VERDANDI_FILE_METHOD_BLUE_HXX


namespace Verdandi
{


    template <class StateErrorVariance, class ObservationOperator,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void ComputeBLUE_matrix(const StateErrorVariance& B,
                            const ObservationOperator& H,
                            const ObservationVector& y,
                            const ObservationErrorVariance& R,
                            StateVector& x);


    template <class StateErrorVariance,
              class ObservationOperator, class CrossedMatrix,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void ComputeBLUE_matrix(const StateErrorVariance& B,
                            const ObservationOperator& H,
                            const CrossedMatrix& cm,
                            const ObservationVector& y,
                            const ObservationErrorVariance& R,
                            StateVector& x);


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BLUE_HXX
#endif
