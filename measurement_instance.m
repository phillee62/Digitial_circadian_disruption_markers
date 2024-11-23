classdef measurement_instance
    properties
        k_observe;
        time_interval;
        hr_phase; 
    end
    methods(Static)
        function obj = initialize()
            obj.k_observe = [];
            obj.time_interval = [];
            obj.hr_phase = [];
        end
        function obj = define(k_observe, time_interval, hr_phase)
            obj.k_observe = k_observe;
            obj.time_interval = time_interval;
            obj.hr_phase = hr_phase;
        end
        function k_observe = return_k_observe(obj)
            k_observe = obj.k_observe;
        end        
        function time_interval = return_time_interval(obj)
            time_interval = obj.time_interval;
        end    
        function hr_phase = return_hr_phase(obj)
            hr_phase = obj.hr_phase;
        end
    end 
end