----------------------------------- GLOBAL -----------------------------------


output_directory = "result/forward-"


----------------------------------- MODEL ------------------------------------


dofile("configuration/quadratic_model.lua")


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      show_iteration = false,
      show_date = true

   },

   output_saver = {

      variable_list = {"state_forecast"},
      file = output_directory .. "%{name}.%{extension}",
      time = "step 100000",
      mode = "binary",
      mode_scalar = "text"

   }

}
