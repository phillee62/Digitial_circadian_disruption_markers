function [hr_phase, phase_std] = import_hr_phase()        

    data = readtable('CRPO_results.csv','Delimiter', ',');
    
    start_day = datenum(table2array(data(1,1)));
    end_day = datenum(table2array(data(end,1)));
    num_sim_days = end_day - start_day + 1;
    
    datestr = [start_day:1:end_day]';
    
    hr_phase = zeros(num_sim_days,1);
    phase_std = zeros(num_sim_days,1);

    for i = 1:num_sim_days
        try
            idx = find(datestr(i) == datenum(table2array(data(:,1))),1);
            phase_est = data(idx,6);
            phase_est = cell2mat(table2array(phase_est));
            
            hr_phase(i) = (datenum(phase_est,'HH:MM')-datenum('00:00','HH:MM'))*24;
            phase_std(i) = table2array(data(idx,7));
        catch 
            hr_phase(i) = NaN;
            phase_std(i) = NaN;
        end
    end
          
end
