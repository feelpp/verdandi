// Copyright (C) 2008-2010, INRIA
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

%include "verdandi.def"

// SWIG interface to Verdandi.
%module verdandi
%{
#define VERDANDI_SPARSE
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "VerdandiHeader.hxx"
#include "VerdandiBase.hxx"

#include "model/QuadraticModel.hxx"
#include "observation_manager/GridToNetworkObservationManager.hxx"
#include "observation_manager/LinearObservationManager.hxx"
#include "method/OptimalInterpolation.hxx"
#include "method/ForwardDriver.hxx"
  %}

%include "std_string.i"
using namespace std;

%import "seldon/seldon.i"

%rename(SeldonError) Seldon::Error;
%rename(OpsError) Ops::Error;

%include "share/Error.hxx"
%include "seldon/share/Errors.hxx"
%include "ops/Error.hxx"

%exception
{
  try
    {
      $action
	}
  catch(Verdandi::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Seldon::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Ops::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "VerdandiHeader.hxx"
%include "VerdandiBase.hxx"

%include "model/QuadraticModel.hxx"
%include "observation_manager/GridToNetworkObservationManager.hxx"
%include "observation_manager/LinearObservationManager.hxx"
%include "method/OptimalInterpolation.hxx"
%include "method/ForwardDriver.hxx"

namespace Verdandi
{
  %template(Model) VSWIG_MODEL;
  %template(GTNObservationManager) VSWIG_GRID_TO_NETWORK_OBSERVATION;
  %template(LObservationManager) VSWIG_LINEAR_OBSERVATION;
  %template(Method) VSWIG_METHOD;
  %template(Method1) VSWIG_METHOD1;
}
