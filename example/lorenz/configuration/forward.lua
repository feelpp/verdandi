----------------------------------- MODEL ------------------------------------


dofile("configuration/lorenz.lua")


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_date = false

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = "result/forward-%{name}.%{extension}",
      mode = "binary",
      mode_scalar = "text"

   }

}
