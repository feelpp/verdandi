----------------------------------- MODEL ------------------------------------


clamped_bar = {

   domain = {

      -- Time step.
      Delta_t = Delta_t_clamped_bar,
      -- Simulation time.
      final_date = final_time_clamped_bar,

      bar_length = 1.,
      Nx = 20

   },

   physics = {

      Young_modulus = 1.,
      mass_density = 1.

   },

   error_statistics = {

      -- Diagonal value of "B".
      background_error_variance = 100.,
      -- Decorrelation length in Balgovind formula.
      background_error_scale = 1.,

   },

   output_saver = {

      variable_list = {"disp_0", "velo_0"},
      file = output_directory .. "/%{name}.bin",
      time = "step " .. Delta_t_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}
