----------------------------------- GLOBAL -----------------------------------


Delta_t_model = 0.0015
Nskip_save = 100

output_directory = "result/"
observation_file = output_directory .. "truth-state_forecast.bin"


----------------------------------- MODEL ------------------------------------


dofile("configuration/quadratic_model.lua")

-- In order to demonstrate the assimilation, errors are introduced in the
-- model.
quadratic_model.definition.initial_state = {0., 0.}
quadratic_model.definition.with_linear_term = false


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "forward-%{name}.%{extension}",
      time = "step " .. Delta_t_model * Nskip_save .. " 1.e-6"

   },

   output = {

      configuration = output_directory .. "forward.lua",
      log = output_directory .. "forward.log"

   }

}


-- Simulation with assimilation using optimal interpolation.
optimal_interpolation = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "vector",

   data_assimilation = {

      analyze_first_step = false,

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "oi-%{name}.%{extension}",
      time = "step " .. Delta_t_model * Nskip_save .. " 1.e-6"

   },

   output = {

      configuration = output_directory .. "oi.lua",
      log = output_directory .. "oi.log"

   }

}


-- Simulation with assimilation using EKF.
extended_kalman_filter = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "matrix",
   -- Computation mode for covariance: "vector" or "matrix".
   covariance_computation = "vector",

   data_assimilation = {

      analyze_first_step = false,

   },

   display = {

      show_iteration = false,
      show_time = true
   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "ekf-%{name}.%{extension}",
      time = "step " .. Delta_t_model * Nskip_save .. " 1.e-6"

   },

   output = {

     configuration = output_directory .. "ekf.lua",
     log = output_directory .. "ekf.log"

  }

}


-- Simulation with assimilation using UKF.
unscented_kalman_filter = {

   data_assimilation = {

      analyze_first_step = false

   },

   sigma_point = {

      -- Choice of sigma-points: "canonical", "star" or "simplex".
      type = "simplex"

   },

   display = {

      show_iteration = false,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "ukf-%{name}.%{extension}",
      time = "step " .. Delta_t_model * Nskip_save .. " 1.e-6"

   },

   output = {

     configuration = output_directory .. "ukf.lua",
     log = output_directory .. "ukf.log"

  }

}
