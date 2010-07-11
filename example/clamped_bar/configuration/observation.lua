-------------------------------- OBSERVATION ---------------------------------


observation = {

   -- Path to the file storing the observations.
   file = observation_file,
   -- How are defined the observations? If the type is "observation, only
   -- observations are stored in the file. If the type is "state", the whole
   -- model state is stored.
   type = "state",
   -- Period with which observations are available.
   Delta_t = Delta_t_clamped_bar * Nskip_save,
   -- Period with which available observations are actually assimilated.
   Nskip = 1,
   -- Duration during which observations are assimilated.
   final_date = final_time_clamped_bar,

   aggregator = {

      type = "step",
      width_left = 0.005,
      width_right = 0.005,

      -- If the value is true, each observation can be used only one time.
      discard_observation = true

   },

   error = {

      -- Variance of observational errors.
      variance = 100.

   },

   operator = {

      -- How is defined the observation operator? If the operator is not
      -- "diagonal", values are read from a file (see entry "File" below).
      definition = "diagonal",
      -- In case of a diagonal operator.
      diagonal_value = 1.,
      -- In case of an operator defined in a file.
      file = "configuration/matrix.bin"

   },

   location = {

      observation_location = {1, 0}

   }

}
