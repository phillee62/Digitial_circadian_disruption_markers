 function LSKF_circadian(dmy_hr, Thr, Tsteps, bin_size)
    
    global model_phase_est model_phase_std fun X_std plotLevelSet
        
    %% Open results.csv file
    fid = fopen('subject.txt');
    results_file = fopen('CRCO_results.csv', 'w');
    fprintf(results_file, 'Date (dd-mmm-yyyy), Weighted Phase Estimate, Weighted Error, Heart Rate Phase Estimate, Heart Rate Phase Error, Model Phase Estimate, Model Phase Error, Raw Model Phase Estimate, Raw Model Phase Error \n');
    
    %% Set up, Import Heart Rate Phase est & std
    % k_noise = 6*10^(-3)*[1 0 0; 0 1 0; 0 0 1]; % System noise 
    k_noise = 6*10^(-2)*[1 0 0; 0 1 0; 0 0 1]; % System noise 

    raw_hr_data = Thr;
    raw_steps_data = Tsteps;
    stp_to_act = mean(Tsteps(:,2));
    [hr_phase_est, hr_phase_std] = import_hr_phase(); % Phase estimation from HR for each day
    num_days = numel(hr_phase_est);
    
    %% Determine whether or not to plot level set
    plotLevelSet = 0;

    %% Initialize variables and initial conditions 
    x_init = [1;0;0.5]; 
%     covariance_init = diag([0.1,0.1,0.1]);  
    covariance_init = diag([0.01,0.01,0.01]);  

    M_init = chol(covariance_init).'; 
    dim = length(covariance_init);
    x_phase = zeros(dim,num_days); % Phase estimation just from model
    M_phase = zeros(dim,dim,num_days);

    model_phase_est = []; 
    model_phase_std = [];

    phase_est = zeros(num_days,1);
    phase_std = zeros(num_days,1);
    
    raw_model_phase_est = zeros(num_days,1);
    raw_model_phase_std = zeros(num_days,1);
    make_fun_raw = [1,0];
    
    fun = griddedInterpolant;
    
    %% Simulate previous day -- i.e, day 0
    try 
        [~, steps_data] = bin_data(raw_hr_data, raw_steps_data, bin_size, 0);
        light = convert_steps_to_light(steps_data(:,2), stp_to_act);
        day0_inst = time_update_instance.define_instance(k_noise, dim, [0:(1/60):24], light, x_init, M_init, [0,0]);
        [x_init, M_init] = sim_day_0(day0_inst);
    catch
        light = [];
        day0_inst = time_update_instance.define_instance(k_noise, dim, [0:(1/60):24], light, x_init, M_init, [0,0]);
        [x_init, M_init] = sim_day_0(day0_inst);   
    end

    raw_model_x_init = x_init;
    raw_model_M_init = M_init;
    
    %% Iterate for each day
    for i = 1:num_days
        
        %% If there is no steps data, make covariance matrix for the next day bigger
        [~, steps_data] = bin_data(raw_hr_data, raw_steps_data, bin_size, i);
        if isempty(steps_data)
            phase_est(i) = NaN;
            phase_std(i) = NaN;
            model_phase_est = [model_phase_est; NaN];
            model_phase_std = [model_phase_std; NaN];
            raw_model_phase_est(i) = NaN;
            raw_model_phase_std(i) = NaN;            
            covar = 2*(M_init*M_init.');
            M_init = chol(covar).';
            
            %% Print results
            fprintf(results_file, '%s, %f, %f, %f, %f, %f, %f, %f, %f\n', datestr(dmy_hr(i+1)),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN);
            
        else
            steps = steps_data(:,2);
            light = convert_steps_to_light(steps, stp_to_act);   
            time_interval = [(i-1)*24*60:1:(i*24*60)]./60;             
            k_observe = hr_phase_std(i)^2; % Measurement noise (from ALSM)
            make_fun = [1,0]; % first entry used to make new interpolation function
                              % second entry used to adjust existing interpolation function after correction
            
            %% Compute raw model output for comparison with LSKF output
            raw_model_time_inst = time_update_instance.define_instance(k_noise,dim,time_interval,light,raw_model_x_init,raw_model_M_init,make_fun_raw);
            [t_phase, raw_xM_phase, raw_xM_init] = raw_model_output(raw_model_time_inst);
            raw_model_phase_est(i) = t_phase;
            raw_model_phase_std(i) = find_phase_std(raw_xM_phase);
            raw_model_x_init = raw_xM_init(:,1);
            raw_model_M_init = raw_xM_init(:,2:4);
            
            %% Modifty HR Phase estimate based on raw model prediction
            hr_phase_est(i) = modify_hr_phase(hr_phase_est(i), t_phase);   
            
            %% Initialize time update & measurement update inst for LSKF
            time_instance = time_update_instance.define_instance(k_noise, dim, time_interval, light, x_init, M_init, make_fun);
            measurement_inst = measurement_instance.define(k_observe,time_interval,hr_phase_est(i));
            
            %% Implement the Level Set Kalman Filter
            [xM_phase, xM_init] = lskf(time_instance, measurement_inst);

            %% Save output, setup for subsequent iteration
            x_phase(:,i) = xM_phase(:,1);
            M_phase(:,:,i) = xM_phase(:,2:4);
            x_init = xM_init(:,1);
            M_init = xM_init(:,2:4);

            %% Compute the phase estimation and error
            phase_std(i) = find_phase_std(xM_phase);
            phase_est(i) = mean(X_std)+1;
            
            %% Print results
            try
                fprintf(results_file, '%s, %s, %f, %s, %f, %s, %f, %s, %f\n', datestr(dmy_hr(i+1)), datestr(phase_est(i)/24, 'HH:MM'), phase_std(i), datestr(hr_phase_est(i)/24, 'HH:MM'), hr_phase_std(i), datestr(model_phase_est(i)/24, 'HH:MM'), model_phase_std(i), datestr(raw_model_phase_est(i)/24, 'HH:MM'), raw_model_phase_std(i));
            catch
                % fprintf(results_file, '%s, %f, %f, %s, %f, %s, %f, %s, %f\n', datestr(dmy_hr(i+1)), NaN, NaN, datestr(hr_phase_est(i)/24, 'HH:MM'), hr_phase_std(i), datestr(model_phase_est(i)/24, 'HH:MM'), model_phase_std(i), datestr(raw_model_phase_est(i)/24, 'HH:MM'), raw_model_phase_std(i));             
                fprintf(results_file, '%s, %f, %f, %s, %f, %s, %f, %s, %f\n', datestr(dmy_hr(i+1)), NaN, NaN, NaN, NaN, datestr(model_phase_est(i)/24, 'HH:MM'), model_phase_std(i), datestr(raw_model_phase_est(i)/24, 'HH:MM'), raw_model_phase_std(i));             
            end

        end
            
    end
end