function [dmy_hr] = ALSM_nonlinear_version(hr_data, steps_data, bin_size, exist_time_zone)

%% Inputs:
% 1. hr_data - a two column array, the first column lists the epoch
% time date in days (e.g., 738157 is 01-Jan-2021), the second column lists
% the heart rate value at that time
% 
% 2. steps_data - a two column array like the heart rate data, but the
% second column is the steps value at the respective time
%
% 3. bin_size - the number of minutes to bin the data (5 is suggested)

%% Outputs:
% A results.csv file that lists the parameter estimates for each day
%
% There are two types of parameter estimate - 1. Mean estiamte of mesor,
% circadian amplitude, circadian phase, and heart rate increase per step
% (HRpS), autocorrelation strength, and 2. Their uncertainty
% estimate 
%
% Moreover, the csv file includes the number of heart rate data points used
% for the estimation, the total step count calculated from steps data
% used for the estimation, the existence of time zone information in
% subject.txt
%
% Here, uncertainty estimate is defined to be the standard deviation of the
% estimated probabilty density function of parameter

%% Global parameters
global hr_avg step_avg N bin_size1

bin_size1 = bin_size;


%% Open the results file

% Open results.csv file and add a header line
results_file = fopen('CRPO_results.csv', 'w');

fprintf(results_file, 'Date (dd-mmm-yyyy), Mean Mesor, Uncertainty of Mesor, Mean Amp., Uncertainty of Amp., Mean Phase, Uncertainty of Phase (hr), Mean HRpS, Uncertainty of HRpS, Mean Autocorr., Uncertainty of Autocorr., Number of Data Points, 2 Day Step Count, Existence of Time Zone Information \n');

%% Process data
hr_dates = hr_data(:,1);
steps_dates = steps_data(:,1);

hr_vals = hr_data(:,2);
steps_vals = steps_data(:,2);

% Save all the days in the interval of data given to loop over
hr_days = floor(hr_dates);
dmy_hr = unique(hr_days, 'rows', 'stable');
dmy_hr = (dmy_hr(1):(dmy_hr(end)+1));
num_days = length(dmy_hr)-2;

%% Remove hr_times and hr_values where steps are zero consecutively for more than 2 hours (i.e., 120 minutes)
[hr_dates, hr_vals] = constant_steps_filter2(hr_dates, hr_vals, steps_dates, steps_vals,120);

%% Main loop

epara=0;

for i = 1:num_days
    
    % Load steps information for that day
    steps_times = steps_dates(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+2));
    steps_values = steps_vals(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+2));
    
    % Load heart rate data for the specific day
    hr_times = hr_dates(hr_dates < dmy_hr(i+2) & hr_dates >= dmy_hr(i));
    hr_values = hr_vals(hr_dates < dmy_hr(i+2) & hr_dates >= dmy_hr(i));
    
    % Convert the times data to minutes
    hr_times = hr_times*24*60;
    steps_times = steps_times*24*60;
    
    step_count = sum(steps_values);
    
    % Ask whether the hr data is empty for the specific 2 days. If so, make
    % the fit values NaN and continue in the loop
    if(isempty(hr_times))
        
        num_points = length(hr_avg(:,1));
        fprintf(results_file, '%s, %f, %f, %f, %f, %s, %f, %f, %f, %f, %f, %i, %i, %f \n', datestr(dmy_hr(i+1)),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, num_points, step_count, exist_time_zone);
        continue;
        
    else
        
        %% THe frist data preprocessing step before the analysis (e.g., Averaging measurements in same bin)
        % This step is the sanme to the preprocessing process described in
        % Bowman et al., 2021.
        raw_steps_data = [steps_times steps_values];
        raw_hr_data = [hr_times hr_values];
        
        % Find all unique binned times for HR data
        [hr_avg, ~, idx] = unique(floor(raw_hr_data(:, 1) / bin_size), 'stable');
        
        % Average HR measurements in same bin
        val = accumarray(idx, raw_hr_data(:, 2), [], @mean);
        
        % Combine unique times and averages
        hr_avg = [hr_avg, val];
        
        try
            
            % Find leftmost and rightmost bin
            left_min = floor(min(raw_hr_data(1, 1), raw_steps_data(1, 1)) / bin_size);
            right_max = floor(max(raw_hr_data(end, 1), raw_steps_data(end, 1)) / bin_size);
            clearvars raw_hr_data;
            
        catch
            
            left_min = floor(raw_hr_data(1, 1) / bin_size);
            right_max = floor(raw_hr_data(end, 1) / bin_size);
            clearvars raw_hr_data;
            
        end
        
        % Total length of interval containing HR data
        period_offset = [left_min, right_max - left_min + 1] * bin_size;
        
        % Number of measurements after averaging same times
        N = size(hr_avg, 1);
        
        % Fill in any gaps in step data with zeros
        step_int = int32(raw_steps_data(:, 1));
        steps_new = zeros(period_offset(2), 2);
        steps_new(:, 1) = period_offset(1) + (0:(period_offset(2) - 1));
        steps_new(step_int - period_offset(1) + 1, 2) = raw_steps_data(:, 2);
        clearvars raw_steps_data;
        
        % Average steps data into bins (including zeros)
        [step_avg, ~, idx] = unique(floor(steps_new(:, 1) / bin_size) , 'stable');
        val = accumarray(idx, steps_new(:, 2), [], @mean);
        step_avg = [step_avg, val];
        
        % We only need to keep the step bins corresponding to HR data
        step_avg = step_avg(hr_avg(:, 1) - step_avg(1, 1) + 1, :);
        
        % If number of measurements after averating same times is smaller than 
        % 20, we do not estimate the parameters, and just return NaN.  
        if(length(hr_avg(:,1))<20)
            
            num_points = length(hr_avg(:,1));
            fprintf(results_file, '%s, %f, %f, %f, %f, %s, %f, %f, %f, %f, %f, %i, %i, %f\n', datestr(dmy_hr(i+1)),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, num_points, step_count, exist_time_zone);
            continue;
            
        end
        
        %% The second data preprocessing step to use ALSM in a suitable manner
        % If a data point has no neighboring data points, it is excluded
        % from the analysis. For instance, there are data points at time t
        % , but there are no data points at time t + binsize (i.e., 5
        % minutes). Then, the points at t are exlcuded from the analysis.
        %
        % Note that this data preprocessing process is required to use the
        % ALSM in a suitable manner.
        
        % arbit_constant: A constant only used for the preprocessing process.
        % Note that you can define the constant to be an arbitrary value
        % which is different to all the wearable measurement values.
        % We recommend using negative values as follows:
        arbit_constant = -10^3; 
       
        time_difference = diff(hr_avg(:,1));
        data_num1 = length(time_difference);
        
        heart_rate1 = arbit_constant*ones(data_num1,1); 
        heart_rate2 = arbit_constant*ones(data_num1,1);
        
        activity1 = arbit_constant*ones(data_num1,1); 
        activity2 = arbit_constant*ones(data_num1,1);
        
        time1 = arbit_constant*ones(data_num1,1); 
        time2 = arbit_constant*ones(data_num1,1);
        
        for ii = 1:data_num1
            
            if abs(time_difference(ii)-1) < 0.01
                
                heart_rate1(ii) = hr_avg(ii,2);
                heart_rate2(ii) = hr_avg(ii+1,2);
                activity1(ii) = step_avg(ii,2);
                activity2(ii) = step_avg(ii+1,2);
                time1(ii) = hr_avg(ii,1) * bin_size;
                time2(ii) = hr_avg(ii+1,1) * bin_size;
                
            end
            
        end
        
        heart_rate1(heart_rate1==arbit_constant)=[]; heart_rate2(heart_rate2==arbit_constant)=[];
        activity1(activity1==arbit_constant)=[]; activity2(activity2==arbit_constant)=[];
        time1(time1==arbit_constant)=[]; time2(time2==arbit_constant)=[];
        
        % If the number of data points after the exclusion descrbied above,
        % is smaller tham 20, we do not estimate the parameters, and just
        % return NaN
        data_num2 = length(heart_rate2);
        
        if data_num2 < 20
            num_points = length(hr_avg(:,1));
            fprintf(results_file, '%s, %f, %f, %f, %f, %s, %f, %f, %f, %f, %f, %i, %i, %f\n', datestr(dmy_hr(i+1)),NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, num_points, step_count, exist_time_zone);
            continue
            
        else
            %% Nonlinear version of Approximation-based least squares method
            
            % we assumed that the free-running period (tau) is 24hr
            % following previous work (Bowman et al., 2021).
            tau = 24; 
            
            % We reformulate and approximate the mathematical model of
            % the heart rate circadian rhythm (Bowman et al., 2021) to
            % convert the model with correlated noise to that with
            % uncorrelated (i.e., Gaussian) noise. Then, we rewrite it in
            % matrix notation.
            xmatrix = zeros(data_num2,5);
            xmatrix(:,1) = cos((2*pi/tau)*mod((time2 / 60), 24));
            xmatrix(:,2) = sin((2*pi/tau)*mod((time2 / 60), 24));
            xmatrix(:,3) = heart_rate1;
            xmatrix(:,4) = activity2;
            xmatrix(:,5) = -activity1;
            
            % creation of initial guess for regression using ALSM
            if epara == 0
                
                % In the first loop, we define the initial guess to be
                % values previously reported in (Bowman et al., 2021).                
                initial_guess = [73, 4*cos((5/4)*pi), 4*sin((5/4)*pi), 0.9, 0.3];
                initial_guess = [initial_guess(1) * (1 - initial_guess(4)), initial_guess(2) * (1 - initial_guess(4)), initial_guess(3) * (1 - initial_guess(4)), initial_guess(4), initial_guess(5)];                
            else
                
                % We use the estimates for day i as the initial guess for
                % day i+1 
                initial_guess = epara';
                
            end
            
            % Reggression model
            % b(1): heart rate mesor * (1- autocorrelation strength)
            % b(2): coefficient of first-order cosine term * (1- autocorrelation strength)
            % b(3): coefficient of first-order sine term * (1- autocorrelation strength)
            % b(4): autocorrelation strength
            % b(5): heart rate increases per one step (HRpS)
            modelfun = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,2) + ...
                b(4)*x(:,3) + b(5)*x(:,4) +b(4)*b(5)*x(:,5);
            
            % Performing a nonlinear regression to estimate the
            % parameters using Levenberg-Marquardt algorithm (Gavin 2020)
            % Here, we simply used the Matlab built-in function to use the
            % algorithm
            estimate = fitnlm(xmatrix, heart_rate2, modelfun, initial_guess);
            
            % Obtaining mean estimates and uncertianty estimates from the
            % regression
            eparainfo = estimate.Coefficients;
            eparainfo = eparainfo{:,:};
            
            epara = eparainfo(1:5,1)'; % Mean estimates
            cov_matrix=estimate.CoefficientCovariance; % Uncertainty estimates
            
            %% Computing the probabilty density function of the parameters
            
            % Here, we approximate the probability density function of
            % circadian parameters (circadian mesor, phase, and amplitude),
            % and HRpS estimate with Monte Carlo Methods
            
            % The number of samples to contruct the density function 
            sample_num = 10^5;
            
            % The distribution of the estimated parameters, b(1), b(2),
            % ..., b(5). Note that the distribution of parameters are
            % assumed to follow Gaussian. This is a fundemental assumption
            % used in regression problems (e.g., Gavin 2020).
            samples=mvnrnd(epara, cov_matrix, sample_num); 
            
            parameter1 = samples(:,1); % Samples of b(1)
            parameter2 = samples(:,2); % Samples of b(2)
            parameter3 = samples(:,3); % Samples of b(3)
            parameter4 = samples(:,4); % Samples of b(4)
            parameter5 = samples(:,5); % Samples of b(5)
            
            % Removing samples having a unrealistic value.
            % 
            % Note that the mesor and HRpS should be greater than 0, and the
            % autocorrelatin strength (alpha) should be between 0 and 1. However,
            % because the support of Gaussian distribution is (-infinite,
            % infinite), their samples can have a unrealistic value.
            % Thus, here, we remove the samples having a unrealistic value.
            
            % Removing mesor sampels having a unrealistic value
            para1=parameter1;
            parameter1(para1 < 0) = [];
            parameter2(para1 < 0) = [];
            parameter3(para1 < 0) = [];
            parameter4(para1 < 0) = [];
            parameter5(para1 < 0) = [];
            
            % Removing alpha samples having a unrealistic value
            para4=parameter4;
            parameter1(para4 < 0) = [];
            parameter2(para4 < 0) = [];
            parameter3(para4 < 0) = [];
            parameter4(para4 < 0) = [];
            parameter5(para4 < 0) = [];
            
            % Note that  if alpha value is very close to 1, unrelastic
            % values of samples of mesor, cosine and sine coefficient can be
            % generated when we obtain them by dividing b(1), b(2), and b(3) by b(4).
            % Thus, here, we use 0.99 instead of 1.
            para41=parameter4;
            parameter1(para41 > 0.99) = [];
            parameter2(para41 > 0.99) = [];
            parameter3(para41 > 0.99) = [];
            parameter4(para41 > 0.99) = [];
            parameter5(para41 > 0.99) = [];
            
            % Removing HRpS sampels having a unrealistic value
            ppara5 = parameter5;
            parameter1(ppara5 < 0) = [];
            parameter2(ppara5 < 0) = [];
            parameter3(ppara5 < 0) = [];
            parameter4(ppara5 < 0) = [];
            parameter5(ppara5 < 0) = [];
            
            samples = [parameter1, parameter2, parameter3, parameter4, parameter5];
                       
            % By dividing b(1), b(2), and b(3) by b(4), we finally obtain
            % sampels of mesor, cosine and sine coefficient.
            
            % The 1st column of tsamples: samples of mesor estimates
            % The 2nd column of tsamples: samples of cosine term estimates
            % The 3rd column of tsampels: samples of sine term estimates
            % The 4th column of tsamples: samples of alpha estimates
            % The 5th column of tsamples: samples of HRpS estimates
            tsamples=[samples(:,1)./(1-samples(:,4)), samples(:,2)./(1-samples(:,4)), ...
                samples(:,3)./(1-samples(:,4)), samples(:,4), samples(:,5)];
            
            % Mesor samples
            mesor_samples = tsamples(:,1);
            mesor_mean = mean(mesor_samples); % Mean estimate of mesor
            mesor_std = std(mesor_samples); % Uncertainty estimate of mesor
            
            % Circadian amplitude samples
            % Here, we compute the amplitude samples using the samples of
            % cosine term and sine term
            amp_samples = sqrt(sum(tsamples(:,2:3).^2, 2));
            amp_mean = mean(amp_samples); % Mean estimate of amplitude
            amp_std = std(amp_samples); % Uncertianty estimate of amplitude
            
            % Alpha samples
            alpha_samples = tsamples(:,4);
            alpha_mean = mean(alpha_samples); % Mean estimate of alpha
            alpha_std = std(alpha_samples); % Uncertainty estimate of alpha
            
            % HRpS samples
            hrps_samples = tsamples(:,5); 
            hrps_mean = mean(hrps_samples); % Mean estimate of HRpS
            hrps_std = std(hrps_samples); % Uncertainty estimate of HRpS
             
            % Phase samples
            % Here, we compute the circadian phase samples using the samples of
            % cosine term and the sine term
            cosvalues = tsamples(:,2) ./ amp_samples;
            sinvalues = tsamples(:,3) ./ amp_samples;
            phase_samples = zeros(length(cosvalues),1);
            
            for kk = 1:length(cosvalues)
                
                sol = angleCalc(sinvalues(kk), cosvalues(kk),'rad');
                sol = sol + pi;
                
                while sol < 0
             
                    sol = sol + 2*pi;
                    
                end
                
                while sol >= 2*pi
                    
                    sol = sol - 2*pi;
                    
                end
                
                phase_samples(kk) = sol;
                
            end
            
            phase_mean = angle(sum(exp(1i*phase_samples)));
            phase_mean = mod(phase_mean, 2*pi);
            phase_est = mod(((24 / (2*pi)) * phase_mean), 24);
            phase_est_str = datestr(phase_est/24, 'HH:MM');
            phase_mean = phase_est_str; % Mean estimate of circadian phase
            
            phase_std = circ_std(phase_samples);
            phase_std = ((24 / (2*pi)) * phase_std); % Uncertianty estimate of circadian phase
            
            %% Write day's results to the results.csv file
            
            num_points = length(hr_avg(:,1));
            fprintf(results_file, '%s, %f, %f, %f, %f, %s, %f, %f, %f, %f, %f, %i, %i, %f\n', datestr(dmy_hr(i+1)), mesor_mean, mesor_std, amp_mean, amp_std, phase_mean, phase_std, hrps_mean, hrps_std, alpha_mean, alpha_std, num_points, step_count, exist_time_zone);
                                                  
        end
        
    end
    
end

end