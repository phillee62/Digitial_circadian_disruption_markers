function [t_phase, raw_xM_phase, raw_xM_init] = raw_model_output(time_instance)
        
    global raw_model

    raw_model = 1;
    [t_phase ,xM_phase, xM_init] = avg_vel_time_update(time_instance);
    raw_xM_phase = xM_phase;
    raw_xM_init = xM_init;
    t_phase = mod(t_phase,24);
end