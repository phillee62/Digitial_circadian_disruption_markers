classdef time_update_instance
    %% Lists properties of time update instance
    properties
        k_noise;
        dim;
        time_interval;   
        light;
        x_mean;
        level_set_span; % This corresponds to the matrix M in the manuscript     
        make_fun;
    end
    
    %% Initialization
    methods(Static)
        function obj = initialize()
            obj.k_noise = [];
            obj.dim = 0;
            obj.time_interval = 0;        
            obj.steps_data = [];
            obj.x_mean = [];
            obj.level_set_span = [];
            obj.make_fun = [1,0];
        end
        function obj = define_instance(k_noise, dim, time_interval, light, x_mean, level_set_span, make_fun)
            obj.k_noise = k_noise;
            obj.dim = dim;
            obj.time_interval = time_interval;
            obj.light = light;
            obj.x_mean = x_mean;
            obj.level_set_span = level_set_span;
            obj.make_fun = make_fun;
        end
        
    %% Class specfific functions
        function obj = modify(obj, tspan, x_mean, level_set_span, make_fun)
            obj.time_interval = tspan;
            obj.x_mean = x_mean;
            obj.level_set_span = level_set_span;
            obj.make_fun = make_fun;
        end
        function [k_noise] = return_k_noise(obj)
            k_noise = obj.k_noise;
        end
        function [dim] = return_dim(obj)
            dim = obj.dim;
        end      
        function [time_step] = return_time_interval(obj)
            time_step = obj.time_interval;
        end
        function [light] = return_light(obj)
            light = obj.light;
        end
        function [x_mean] = return_x_mean(obj)
            x_mean = obj.x_mean;
        end
        function [level_set_span]= return_level_set_span(obj)
            level_set_span = obj.level_set_span;
        end
    end
end