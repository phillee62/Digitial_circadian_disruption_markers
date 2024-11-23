function [t_phase, xM_phase, xM_init] = avg_vel_time_update(time_instance)
    
    global fun plotLevelSet

    k_noise = time_instance.k_noise; 
    light = time_instance.light;
    tspan = time_instance.time_interval;
    x_mean = time_instance.x_mean;
    M = time_instance.level_set_span;     
    x_and_M = [x_mean, M];   
    new_fun = time_instance.make_fun(1);
    modify_fun = time_instance.make_fun(2);
    
    try 
        %% Propagate the system forward
        [t,y] = ode45(@(t,y)dxMdt(t,y,x_and_M,k_noise,light),tspan,x_and_M);
        minx_idx = find(y(:,1) == min(y(:,1)));
        circ_phase = y(minx_idx,:);
        t_phase = t(minx_idx);
        xM_phase = reshape(circ_phase,size(x_and_M));
        xM_init = reshape(y(numel(t),:),size(x_and_M));
        
        %% Create or modify interpolation function                    
        if new_fun && ~modify_fun
            fun = interpolate_data(t,y);  
        elseif ~new_fun && modify_fun
            modify_measurement_fun(t,y);
        end
        
    catch
        t_phase = tspan(1);
        xM_phase = x_and_M;
        xM_init = x_and_M;
    end

end