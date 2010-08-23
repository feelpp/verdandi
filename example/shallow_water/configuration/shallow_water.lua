----------------------------------- MODEL ------------------------------------


shallow_water = {

   domain = {

      -- Time step.
      Delta_t = Delta_t_shallow_water,
      -- Simulation time.
      final_time = final_time_shallow_water,

      x_min = 0.,
      y_min = 0.,

      Delta_x = 1.,
      Delta_y = 1.,

      Nx = 100,
      Ny = 1

   },

   initial_condition = {

      -- Value that may be taken in the center or on the left.
      value = 1.05,
      -- Take 'value' in the center?
      center = true,
      -- Take 'value' on the left?
      left = false

   },

   model_error = {

      -- Standard deviation of the white Gaussian noise on *boundary*
      -- conditions.
      standard_deviation_bc = 0.,
      -- Standard deviation of the white Gaussian noise on *initial*
      -- conditions.
      standard_deviation_ic = 0.,

      -- Newran seed directory (with "/" at the end), or "current_time" for
      -- random seeds generated with current time, or a given seed number (in
      -- ]0, 1[).
      random_seed = 0.5

   },

   error_statistics = {

      -- Diagonal value of "B".
      state_error_variance = 100.,
      -- Decorrelation length in Balgovind formula.
      state_error_scale = 1.,

      -- Not used in our case.
      model_error_variance = 1.,
      -- Not used in our case.
      model_error_scale = 1.

   },

   boundary_condition = {

      -- Available options:
      --   free
      --   wall
      --   flow [value]
      --   height [value]
      left = "flow 0.1",
      right = "free",
      bottom = "wall",
      top = "wall"

   },

   data_assimilation = {

      with_positivity_requirement = true

   },

   output_saver = {

      variable_list = {"u", "v", "h"},
      file = output_directory .. "/%{name}.bin",
      time = "step " .. Delta_t_shallow_water * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   }

}
