----------------------------------- GLOBAL -----------------------------------


Delta_t_shallow_water = 0.03
final_time_shallow_water = 1500 * Delta_t_shallow_water
-- Saving period.
Nskip_save = 10

output_directory = "result/river-truth/"
output_mode = "binary"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/shallow_water.lua")


----------------------------------- METHOD -----------------------------------


forward = {

   display = {

      show_iteration = false,
      show_date = false

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "/forward-%{name}.%{extension}",
      time = "step " .. Delta_t_shallow_water * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}
