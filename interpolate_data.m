function fun = interpolate_data(t,y)
    
    t_start = t(1);
    t = t-t_start;
    angle_sample = zeros(numel(t),1);

    for i=1:numel(t)
       angle_sample(i) = convert_to_angle(y(i,1:2));
    end

    sorted_ang = sort(angle_sample);
    sorted_t = zeros(numel(t),1);
    
    for j = 1:numel(t)
        sorted_t(find(angle_sample(j) == sorted_ang,1)) = t(j);
    end
    
    fun = griddedInterpolant(sorted_ang, sorted_t);
end
