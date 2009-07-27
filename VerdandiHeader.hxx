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

#include "seldon/SeldonHeader.hxx"
#include "talos/TalosHeader.hxx"


namespace Verdandi
{


    using namespace std;

    using namespace Seldon;
    using namespace Talos;

    using Seldon::to_num;
    using Seldon::to_str;


} // namespace Verdandi.


#include "share/Error.hxx"
#include "share/UsefulFunction.hxx"


#define VERDANDI_FILE_VERDANDIHEADER_HXX
#endif
