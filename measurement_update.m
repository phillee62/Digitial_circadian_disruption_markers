function [x_corrected, M_corrected] = measurement_update(measurement_inst, xM_phase)
    
%   Input Arguments:
%   xM_phase: average and a square root of predicted covariance matrix of class mean_covariance_sqrt_cls
%   measurement: heart rate phase measurement
%   measurement_fun: (non-stochastic) measurement function z = f(x)
%   measurement_noise: Covariance matrix of measurement noise (not a square root)

%   Output Argument:
%   x_corrected: corrected state vector at CBT min
%   M_corrected: corrected covariance matrix at CBT min
    
    measurement = measurement_inst.hr_phase;
    measurement_fun = @measurement_function;
    measurement_noise = measurement_inst.k_observe;

   [x_corrected, M_corrected] = square_root_cubature_measurement_update(xM_phase, measurement, measurement_fun, measurement_noise);
 
end