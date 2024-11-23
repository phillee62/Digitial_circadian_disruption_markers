function phase = measurement_function(x_mean)
    
    global fun

    angle = convert_to_angle(x_mean);
    phase = fun(angle)-1;
    
end