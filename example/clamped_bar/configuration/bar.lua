----------------------------------- GLOBAL -----------------------------------


Delta_t_clamped_bar = 0.01
final_time_clamped_bar = 10.0
-- Saving period.
Nskip_save = 10

observation_file = "result/bar-truth/forward-state_forecast.bin"

output_directory = "result/bar/"
output_mode = "binary"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/clamped_bar.lua")


-- The configuration of the model is modified.
clamped_bar.physics.mass_density = 1.


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")


----------------------------------- METHOD -----------------------------------


-- Simulation with assimilation using optimal interpolation.
optimal_interpolation = {

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "matrix",

   data_assimilation = {

      analyze_first_step = true,

   },

   display = {

      show_iteration = false,
      show_date = false
   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "/oi-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_date = false

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "/forward-%{name}.%{extension}",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}
