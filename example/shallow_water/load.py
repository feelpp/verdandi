# This script loads the results of the model. It defines:
#   Nx: number of points along x
#   Ny: number of points along y
#   h: water height   [dimensions: Nt, Nx, Ny]
#   u: velocity along x   [Nt, Nx, Ny]
#   v: velocity along y   [Nt, Nx, Ny]
# If available:
#   h_a: analyzed water height   [Nt, Nx, Ny]
#   Ne: number of members in the ensemble
#   he: ensemble water heights   [Ne, Nt, Nx, Ny]
#   he_a: analyzed ensemble water heights   [Ne, Nt, Nx, Ny]
#   Np: number of points where the cost function is evaluated with program
#       "eval_cost".
#   eval_cost_point: point where the cost function is evaluated by program
#                    "eval_cost"   [Np]
#   eval_cost: cost function from program "cost"   [Np]
#   eval_cost_bg: background part of the cost function, program "cost"   [Np]
#   eval_cost_obs: observation part of cost function, program "cost"   [Np]
#   Niter: number of iterations in 4D-Var
#   cost: cost function from 4D-Var   [Niter]
#   cost_gradient: gradient of cost function from 4D-Var   [Niter, Nx, Ny]
#   cost_bg: background part of the cost function, 4D-Var   [Niter]
#   cost_obs: observation part of cost function, 4D-Var   [Niter]

import os, sys, glob
from atmopy.talos import *
from numpy import *

if len(sys.argv) != 2:
    raise Exception, "One (and only one) directory needed!"
result_directory = sys.argv[1]
if not os.path.isdir(result_directory):
    raise Exception, "\"" + result_directory + "\" is not a directory!"

# Configuration.
configuration_file = glob.glob(result_directory + "/*.cfg")
if len(configuration_file) == 0:
    raise Exception, "Unable to find the configuration file (.cfg) in \"" \
        + result_directory + "\"."
if len(configuration_file) > 1:
    raise Exception, "Found several configuration files (.cfg) in \"" \
        + result_directory + "\". Please remove the unrelated configuration" \
        " files."
configuration_file = configuration_file[0]

# Dimensions of the domain.
configuration = Config(configuration_file)
Nx = configuration.Nx
Ny = configuration.Ny

# Reading the main data, if available.
if os.path.isfile(os.path.join(result_directory, "h.bin")):
    h = fromfile(os.path.join(result_directory, "h.bin"), dtype = "Float64")
    h.shape = (h.shape[0] / (Nx * Ny), Nx, Ny)
    Nt = h.shape[0]
    u = fromfile(os.path.join(result_directory, "u.bin"), dtype = "Float64")
    u.shape = (Nt, Nx, Ny)
    v = fromfile(os.path.join(result_directory, "v.bin"), dtype = "Float64")
    v.shape = (Nt, Nx, Ny)
# Analysis also, if any.
if os.path.isfile(os.path.join(result_directory, "h_a.bin")):
    h_a = fromfile(os.path.join(result_directory, "h_a.bin"),
                   dtype = "Float64")
    h_a.shape = (h_a.shape[0] / (Nx * Ny), Nx, Ny)
# Ensemble forecasts also, if any.
if os.path.isfile(os.path.join(result_directory, "h-000.bin")):
    he = []
    file_list = glob.glob(os.path.join(result_directory, "h-*.bin"))
    for filename in sort(file_list):
        tmp = fromfile(filename, dtype = "Float64")
        tmp.shape = (tmp.shape[0] / (Nx * Ny), Nx, Ny)
        he.append(tmp)
    he = array(he)
    Ne = he.shape[0]
# Ensemble analyses also, if any.
if os.path.isfile(os.path.join(result_directory, "h_a-000.bin")):
    he_a = []
    file_list = glob.glob(os.path.join(result_directory, "h_a-*.bin"))
    for filename in sort(file_list):
        tmp = fromfile(filename, dtype = "Float64")
        tmp.shape = (tmp.shape[0] / (Nx * Ny), Nx, Ny)
        he_a.append(tmp)
    he_a = array(he_a)
    Ne = he_a.shape[0]
# Cost function, if available (program "eval_cost").
if os.path.isfile(os.path.join(result_directory, "eval_cost_point.bin")):
    eval_cost = fromfile(os.path.join(result_directory, "eval_cost.bin"),
                         dtype = "Float64")
    eval_cost_bg = fromfile(os.path.join(result_directory,
                                         "eval_cost_bg.bin"),
                            dtype = "Float64")
    eval_cost_obs = fromfile(os.path.join(result_directory,
                                          "eval_cost_obs.bin"),
                             dtype = "Float64")
    eval_cost_point = fromfile(os.path.join(result_directory,
                                            "eval_cost_point.bin"),
                               dtype = "Float64")
    Np = eval_cost_point.shape[0]
# 4D-Var outputs.
if os.path.isfile(os.path.join(result_directory, "cost.bin")):
    cost = fromfile(os.path.join(result_directory, "cost.bin"),
                    dtype = "Float64")
    cost_bg = fromfile(os.path.join(result_directory, "cost_bg.bin"),
                       dtype = "Float64")
    cost_obs = fromfile(os.path.join(result_directory, "cost_obs.bin"),
                        dtype = "Float64")
    cost_gradient = fromfile(os.path.join(result_directory,
                                          "cost_gradient.bin"),
                             dtype = "Float64")
    Niter = cost.shape[0]
    cost_gradient.shape = (Niter, Nx, Ny)
    ic = fromfile(os.path.join(result_directory, "ic.bin"), dtype = "Float64")
    ic.shape = (Niter + 1, Nx, Ny)

# Convenient functions.
from atmopy.stat import spatial_distribution
from atmopy.stat import time_evolution
ensemble_measure = spatial_distribution
def rmse(a, b):
    return sqrt(((a - b)**2).mean())
def rms(a):
    return sqrt((a**2).mean())
