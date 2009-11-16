// Copyright (C) 2008-2009 INRIA
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


#ifndef VERDANDI_FILE_VERDANDIHEADER_HXX


//! The namespace of the data assimilation library Verdandi.
namespace Verdandi
{


} // namespace Verdandi.


#include <iostream>
#include <fstream>
#include <map>
#include <list>

#include "seldon/SeldonHeader.hxx"
#include "getpot/GetPot.hpp"


namespace Verdandi
{


    using namespace std;

    using namespace Seldon;

    using Seldon::to_num;
    using Seldon::to_str;


} // namespace Verdandi.


#include "share/Error.hxx"
#include "share/UsefulFunction.hxx"
#include "share/MessageHandler.hxx"

#include "output_saver/OutputSaver.hxx"

#ifdef VERDANDI_SPARSE
#define VERDANDI_TANGENT_OPERATOR_SPARSE
#define VERDANDI_OBSERVATION_ERROR_SPARSE
#define VERDANDI_BACKGROUND_ERROR_SPARSE
#endif

#ifdef VERDANDI_DENSE
#define VERDANDI_TANGENT_OPERATOR_DENSE
#define VERDANDI_BACKGROUND_ERROR_DENSE
#endif

#define VERDANDI_FILE_VERDANDIHEADER_HXX
#endif
