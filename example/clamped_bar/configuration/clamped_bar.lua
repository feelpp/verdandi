----------------------------------- MODEL ------------------------------------


clamped_bar = {

   domain = {

      -- Time step.
      Delta_t = Delta_t_clamped_bar,
      -- Simulation time.
      final_time = final_time_clamped_bar,

      bar_length = 1.,
      Nx = 20

   },

   physics = {

      Young_modulus = 1.,
      mass_density = 1.

   },

   state_error = {

      -- Is the state error variance a scaled identity matrix?
      scaled_identity = true,
      -- If so, put the diagonal value:
      diagonal_value = 100.,
      -- Otherwise, the operator value (file name or table):
      value = {}

   },

   output_saver = {

      variable_list = {"disp_0", "velo_0"},
      file = output_directory .. "/%{name}.bin",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}
