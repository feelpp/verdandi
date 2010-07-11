----------------------------------- MODEL ------------------------------------


lorenz = {

   parameter = {

      Prandtl = 10.,
      Rayleigh = 28.,
      b = 2.6666667

   },

   initial_condition = {

      X = 0.1,
      Y = 0.1,
      Z = 0.1

   },

   time = {

      Delta_t = .001,
      initial_date = 0.,
      final_date = 100.

   },

   output_saver = {

      variable_list = {"X", "Y", "Z"},
      file = "result/forward-%{name}.bin",
      mode_scalar = "binary"

   }

}
