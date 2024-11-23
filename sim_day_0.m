function [x_init, M_init] = sim_day_0(day0_inst)
    
    light = day0_inst.light;
    
    if ~isempty(light)
        [~, ~, xM_init] = avg_vel_time_update(day0_inst);
        x_init = xM_init(:,1);
        M_init = xM_init(:,2:4);
    else
        x_init = day0_inst.x_mean;
        M_init = day0_inst.level_set_span;
    end
end