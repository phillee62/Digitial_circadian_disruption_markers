sleep_data = readtable('Sleep_results.csv','Delimiter', ',');
data = readtable('CRCO_results.csv','Delimiter', ',');

start_day = datenum(table2array(data(1,1)));
end_day = datenum(table2array(data(end,1)));
num_sim_days = end_day - start_day + 1;

datestring = [start_day:1:end_day]';

CRPO_phase = zeros(num_sim_days,1);
CRCO_phase = zeros(num_sim_days,1);

CRCO_phase_before = zeros(num_sim_days,1);
CRCO_uncertainty_before = zeros(num_sim_days,1);

sleep_midpt = zeros(num_sim_days,1);

for i = 1:num_sim_days
    try
        %% Import estimated CRCO and CRPO phase from "CRCO_results.csv"
        % to calculate CRCO-sleep and CRPO-sleep misalignment

        idx = find(datestring(i) == datenum(table2array(data(:,1))),1);
        CRCO_est = data(idx,2);
        CRCO_est = cell2mat(table2array(CRCO_est));
        CRPO_est = data(idx,4);
        CRPO_est = cell2mat(table2array(CRPO_est));
        
        CRCO_phase(i) = (datenum(CRCO_est,'HH:MM')-datenum('00:00','HH:MM'))*24;
        CRPO_phase(i) = (datenum(CRPO_est,'HH:MM')-datenum('00:00','HH:MM'))*24;

        %% Import (before correction in Kalman Filtering) CRCO and CRPO uncertainties from "CRCO_results.csv"
        % to calculate Internal misalignment

        CRCO_before_correction = data(idx,6);
        CRCO_before_correction = cell2mat(table2array(CRCO_before_correction));
        CRCO_before_correction_uncertainty = data(idx,7);
        CRCO_before_correction_uncertainty = table2array(CRCO_before_correction_uncertainty);
        
        CRCO_phase_before(i) = (datenum(CRCO_before_correction,'HH:MM')-datenum('00:00','HH:MM'))*24;
        CRCO_uncertainty_before(i) = CRCO_before_correction_uncertainty;


        %% Import estimated sleep midpoint from "Sleep_results.csv"
        % to calculate internal misalignment

        idx_sleep = find(datestring(i) == datenum(table2array(sleep_data(:,1))),1);
        midpt = sleep_data(idx_sleep,2);
        midpt = table2array(midpt);
        sleep_midpt(i) = midpt;


    catch 
        CRCO_phase(i) = NaN;
        CRPO_phase(i) = NaN;
        CRCO_phase_before(i) = NaN;
        CRCO_uncertainty_before(i) = NaN;
        sleep_midpt(i) = NaN;
    end
end

%% Compute each disruption marker
phi_ref = -1;

CRCO_misalign = (CRCO_phase + phi_ref) - sleep_midpt;
CRCO_misalign(CRCO_misalign > 12) = CRCO_misalign(CRCO_misalign > 12) - 24;
CRCO_misalign(CRCO_misalign < -12) = CRCO_misalign(CRCO_misalign < -12) + 24;
CRCO_misalign = abs(CRCO_misalign);

CRPO_misalign = CRPO_phase - sleep_midpt;
CRPO_misalign(CRPO_misalign > 12) = CRPO_misalign(CRPO_misalign > 12) - 24;
CRPO_misalign(CRPO_misalign < -12) = CRPO_misalign(CRPO_misalign < -12) + 24;
CRPO_misalign = abs(CRPO_misalign);

Internal_misalign = (CRCO_phase_before + phi_ref) - CRPO_phase;
Internal_misalign(Internal_misalign > 12) = Internal_misalign(Internal_misalign > 12) - 24;
Internal_misalign(Internal_misalign < -12) = Internal_misalign(Internal_misalign < -12) + 24;
Internal_misalign = abs(Internal_misalign)./CRCO_uncertainty_before;

%% Save output in .csv file.
results_file = fopen('Digital_circadian_disruption_markers.csv', 'w');
fprintf(results_file, 'Date (dd-mmm-yyyy), CRCO-sleep_misalignment, CRPO-sleep_misalignment, Internal_misalignment \n');

date_column = string(table2array(data(:,1)));
for j = 1:num_sim_days
    fprintf(results_file, '%s, %f, %f, %f\n', date_column(j), CRCO_misalign(j), CRPO_misalign(j), Internal_misalign(j));
end

