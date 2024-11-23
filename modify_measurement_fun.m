function modify_measurement_fun(t,y)
    
    global fun
    
    t = mod(t,24);
    t_prev = fun.Values;
    ang_prev = cell2mat(fun.GridVectors);
    ang = zeros(length(y),1);
    
    for i = 1:length(ang)
        ang(i) = convert_to_angle(y(i,1:2));
    end
    
    t_phase = mod(t(1),24);
    ang_phase = convert_to_angle(y(1,1:2));
    
    idx = find(t_prev < t_phase);
    
    t_new = [t_prev(idx); t];
    ang_new = [ang_prev(idx); ang];
    
    sorted_ang = sort(ang_new);
    sorted_t = zeros(numel(t_new),1);
    
    for j = 1:numel(t_new)
        sorted_t(find(ang_new(j) == sorted_ang,1)) = t_new(j);
    end
    
    fun = griddedInterpolant(sorted_ang, sorted_t);
    

end