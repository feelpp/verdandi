----------------------------------- GLOBAL -----------------------------------


Delta_t_shallow_water = 0.03
final_time_shallow_water = 1500 * Delta_t_shallow_water
-- Saving period.
Nskip_save = 10

observation_file = "result/river-truth/forward-state_forecast.bin"

output_directory = "result/river/"
output_mode = "binary"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/shallow_water.lua")

-- The configuration of the model is modified.
shallow_water.initial_condition.value = 1.
shallow_water.error.standard_deviation_bc = 0.1


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")


----------------------------------- METHOD -----------------------------------


-- Simulation with assimilation using optimal interpolation.
optimal_interpolation = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "vector",

   data_assimilation = {

      analyze_first_step = true,

   },

   display = {

      show_iteration = false,
      show_time = false
   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "oi-%{name}.%{extension}",
      time = "step " .. Delta_t_shallow_water * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "oi.lua",
      log = output_directory .. "oi.log"

   }

}


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_time = false

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "forward-%{name}.%{extension}",
      time = "step " .. Delta_t_shallow_water * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "forward.lua",
      log = output_directory .. "forward.log"

   }

}
