function run_sleep_analysis(Tsleep, Tsteps, isEmptySleep, isEmptySteps, exist_time_zone)

    if(~isEmptySleep) && (~isEmptySteps)
        min_sleep_duration = 2; % minimum sleep duration
        wbstime = 20/60; % neglibible wakeness time betwen sleeps
        steps_times = Tsteps(:,1); % steps times 
        steps_values = Tsteps(:,2); % steps counts
        sleep_step_filter(Tsleep, min_sleep_duration, wbstime, steps_times, steps_values, exist_time_zone);

    elseif (~isEmptySleep) 
        min_sleep_duration = 2; % minimum sleep duration
        wbstime = 20/60; % neglibible wakeness time betwen sleeps
        sleep_filter(Tsleep, min_sleep_duration, wbstime, exist_time_zone);

    elseif (~isEmptySteps)
        min_sleep_duration = 2; % minimum sleep duration
        wbstime = 20/60; % neglibible wakeness time betwen sleeps
        steps_times = Tsteps(:,1); % steps times 
        steps_values = Tsteps(:,2); % steps counts
        step_filter(min_sleep_duration, wbstime, steps_times, steps_values, exist_time_zone);
    end
    
end