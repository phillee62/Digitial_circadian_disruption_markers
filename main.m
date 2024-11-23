%% This is the main execution file.
clc; clear; close all;

%% Load the subject file to get the timezone information, if available
fid = fopen('Sample_data/subject.txt');
tline = fgetl(fid);
sub_file = cell(0,1);
while ischar(tline)
    sub_file{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);
timezone = sub_file{2};

ref_date = datenum(1970,1,1);

%% Average werable data in 5-minute bins are analyzed following previous work (Bowman et al., 2021, Kim et al., 2023)
bin_size = 5;

% Existance of time zone information in subject.txt
% exist_time_zone = 1 means that the information exist in subject. text, 
% and exist_time zone = 0 represents that it does not exist 

try
    if(isempty(timezone))
	t_offset=0;
    exist_time_zone = 0;
    else
    	t_offset = tzoffset(datetime('today', 'TimeZone', timezone));
    	t_offset = datenum(t_offset);
        exist_time_zone = 1;
    end
catch
    t_offset = 0;
    exist_time_zone = 0;
end

%% Load wearable heart rate data
try
    Thr = readtable('Sample_data/heart_rate.csv','Delimiter', ',');
    Thr = table2array(Thr(:,1:2));
    Thr(:,1) = Thr(:,1)/(24*60*60) + ref_date + t_offset;
catch
    error('No heart rate file in directory.')
end

%% Load wearable steps data
try
    Tsteps = readtable('Sample_data/steps.csv','Delimiter', ',');
    Tsteps = table2array(Tsteps(:,1:2));
    Tsteps(:,1) = Tsteps(:,1)/(24*60*60) + ref_date + t_offset;
    isEmptySteps = 0;
catch
    Tsteps = [];
    isEmptySteps = 1;
end

% Load weaerable sleep data describing whether a subject is awake for each
% time
try
    Tsleep = readtable('Sample_data/sleep.csv','Delimiter', ',');
    Tsleep = table2array(Tsleep(:,1:2));
    Tsleep(:,1) = Tsleep(:,1)/(24*60*60) + ref_date + t_offset;
    isEmptySleep = 0;
catch
    Tsleep = [];
    isEmptySleep = 1;
end
if(~isEmptySleep)
    Thr = remove_sleep(Thr, Tsleep);
end
if(~isEmptySteps)
    %% Estimate sleep onset, offset, and sleep offset based on available sleep activity data
    run_sleep_analysis(Tsleep, Tsteps, isEmptySleep, isEmptySteps, exist_time_zone);

    %% Circadian Rhythms in the Peripheral Oscillator (CRPO; Kim et al., J. R. Soc. Interface, 2023)
    dmy_hr = ALSM_nonlinear_version(Thr, Tsteps, bin_size, exist_time_zone);
    
    %% Analyze Circadian Rhythms in the Central Oscillator (CRCO; Kim et al., SIAM J. Appl. Math, 2024) 
    LSKF_circadian(dmy_hr, Thr, Tsteps, bin_size);

    %% Compute the digital markers of circadian disruptions
end