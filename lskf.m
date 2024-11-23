function [xM_phase, xM_init] = lskf(time_instance, measurement_inst)

    global model_phase_est model_phase_std
        
    k_observe = measurement_instance.return_k_observe(measurement_inst);
    hr_phase = measurement_instance.return_hr_phase(measurement_inst);
    t_start = time_instance.time_interval(1);
    t_end = time_instance.time_interval(end);

    if isnan(k_observe) || isnan(hr_phase)
        %% Time update (No measurement update - i.e. no HR, yes Steps data)
        [t_phase, xM_phase, xM_init] = avg_vel_time_update(time_instance);
        phase = t_phase - t_start;
        model_phase_est = [model_phase_est; phase];
        model_phase_std = [model_phase_std; find_phase_std(xM_phase)];
        
    else
        %% Time update (Before phase correction)
        [t_phase, xM_phase, ~] = avg_vel_time_update(time_instance);
        
        %% Save predicted phase estimate and error          
        model_phase_std = [model_phase_std; find_phase_std(xM_phase)];
        phase = mod(t_phase, 24);
        model_phase_est = [model_phase_est; phase];
        
        %% Measurement update step    
        [x_corrected, M_corrected] = measurement_update(measurement_inst, xM_phase);    
        xM_phase = [x_corrected, M_corrected];
        
        %% Time update (After phase correction)
        t_pha = measurement_function(x_corrected)+1;        
        phase_corrected = t_pha + t_start;        
        t_int = [phase_corrected:(1/60):t_end];
        make_fun = [0,1];
        time_instance = time_update_instance.modify(time_instance, t_int, x_corrected, M_corrected, make_fun);
        [~, ~, xM_init] = avg_vel_time_update(time_instance);

    end
end
