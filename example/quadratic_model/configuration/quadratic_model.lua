----------------------------------- MODEL ------------------------------------


quadratic_model = {

   definition = {

      initial_state = {0.3,
                       0.3},

      with_quadratic_term = true,
      with_linear_term = true,
      with_constant_term = true,

      quadratic_term = {1., 0.,
                        0., 1.,
                        1., 0.,
                        0., 0.},
      linear_term = {-1., 1.,
                     0., -1.},
      constant = {0.,
                  0.},

      Delta_t = .0015,
      initial_date = 0.,
      final_date = 2.5

   },

   output_saver = {

      variable_list = {"state", "Q", "L", "b"},
      file = output_directory .. "%{name}.bin",
      time = "step 0.015 1.e-6"

   }

}
