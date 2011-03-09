----------------------------------- MODEL ------------------------------------


parametric_clamped_bar = {

    -- Components of the state vector.
    state = {"displacement", "velocity", "theta_force"},
    -- Reduced state parameters.
    reduced_state = {"theta_force"},

    domain = {

        -- Time step.
        Delta_t = Delta_t_parametric_clamped_bar,
        -- Simulation time.
        final_time = final_time_parametric_clamped_bar,

        bar_length = 1.,
        Nx = 100

    },

    physics = {

        Young_modulus = 1.,
        mass_density = 1.,
        theta_force = {1., 1., 1., 1.},
        theta_stiffness = {1., 1.},
        theta_mass = {1., 1., 1., 1., 1.},
        theta_damp = {1.},
        alpha = 1.,
        beta = 1.
    },

    error_statistics = {

        -- Diagonal value of "B".
        state_error_variance = 100.,
        -- Decorrelation length in Balgovind formula.
        state_error_scale = 1.,

    },

    output_saver = {

        variable_list = {"disp_0", "velo_0"},
        file = output_directory .. "/%{name}.bin",
        time = "step " .. Delta_t_parametric_clamped_bar * Nskip_save .. " 1.e-6",
        mode = output_mode,
        mode_scalar = output_mode_scalar

   }

}
