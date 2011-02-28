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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_VERDANDIHEADER_HXX


//! The namespace of the data assimilation library Verdandi.
namespace Verdandi
{


} // namespace Verdandi.


#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <cmath>


// Convenient macros to catch exceptions.
#ifndef TRY
#define TRY try {
#endif
#ifndef END
#define END                                                     \
  }                                                             \
    catch(Verdandi::Error& Err)                                 \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(Seldon::Error& Err)                                   \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(::Ops::Error& Err)                            \
      {                                                         \
        Err.CoutWhat();                                         \
        return 1;                                               \
      }                                                         \
    catch(std::exception& Err)                                  \
      {                                                         \
        cout << "C++ exception: " << Err.what() << endl;        \
        return 1;                                               \
      }                                                         \
    catch(std::string& str)                                     \
      {                                                         \
        cout << str << endl;                                    \
        return 1;                                               \
      }                                                         \
    catch(const char* str)                                      \
      {                                                         \
        cout << str << endl;                                    \
        return 1;                                               \
      }                                                         \
    catch(...)                                                  \
      {                                                         \
        cout << "Unknown exception..." << endl;                 \
        return 1;                                               \
      }
#endif

#ifdef VERDANDI_WITH_ABORT
#define SELDON_WITH_ABORT
#endif

#ifndef VERDANDI_WITH_ABORT
#define OPS_WITH_EXCEPTION
#endif

#include "seldon/SeldonHeader.hxx"
#include "seldon/vector/Vector2.hxx"
#include "seldon/vector/Vector3.hxx"


namespace Verdandi
{


    using namespace std;

    using namespace Seldon;

    using Seldon::to_num;
    using Seldon::to_str;


} // namespace Verdandi.


#include "share/Ops.hxx"
#include "share/Logger.hxx"
#include "share/Error.hxx"
#include "share/UsefulFunction.hxx"
#include "share/MessageHandler.hxx"
#include "share/VerdandiBase.hxx"

#include "output_saver/OutputSaver.hxx"

#ifdef VERDANDI_SPARSE
#define VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
#define VERDANDI_OBSERVATION_ERROR_SPARSE
#define VERDANDI_STATE_ERROR_SPARSE
#endif

#ifdef VERDANDI_DENSE
#define VERDANDI_TANGENT_LINEAR_OPERATOR_DENSE
#define VERDANDI_STATE_ERROR_DENSE
#endif

#define VERDANDI_FILE_VERDANDIHEADER_HXX
#endif
