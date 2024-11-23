function sleep_step_filter(Tsleep, min_sleep_duration, wbstime, steps_times, steps_values, exist_time_zone)

%%% Function that calculates the midpoint of sleep.

% Open results.csv file and add a header line
results_file = fopen('Sleep_results.csv', 'w');
fprintf(results_file, 'Date (dd-mmm-yyyy), sleep_midpoint, sleep_onset, sleep_offset, Existence_of_Time_Zone_Information \n');

%% Estimation of sleep phase from an activity file

% If steps are zero csecutively for more than min_sleep_duration, it is
% consdiered as a sleep phase

s_length = length(steps_times);

is_zero = 0;

starts = [];
ends = [];

% Find all segments of zeros

for i = 1:s_length
    
    if(~is_zero)
        
        if (i < s_length) && (~steps_values(i))
            
            starts = [starts; steps_times(i)];
            
            is_zero = 1;
            
        end
        
    else
        
        if(~steps_values(i))
            
            if(i == s_length)
                ends = [ends; steps_times(end)];
            else
                
                continue;
                
            end
            
        else
            
            ends = [ends; steps_times(i-1)];
            
            is_zero = 0;
            
        end
        
    end
    
end

segment_lengths = (ends-starts)*24;
v = find(segment_lengths>min_sleep_duration);
%sleep_segments = [starts(v) ends(v)];

filter_start_times = starts(v);
filter_end_times = ends(v);

filter_dmy_sleep = floor(filter_start_times);

index = 1;
filter_data = {};

for j = 1:length(filter_dmy_sleep)
    
    if j == index
        
        start_time_collection = [];
        end_time_collection = [];
        sleep_collection =[];
        
        start_time_collection = [start_time_collection; filter_start_times(j,1)];
        end_time_collection = [end_time_collection; filter_end_times(j,1)];
        sleep_collection = [sleep_collection; filter_end_times(j,1) - filter_start_times(j,1)];
        
        if j < length(filter_dmy_sleep)
            
            wcount = 0;
            
            while index < length(filter_dmy_sleep) && filter_dmy_sleep(j) == filter_dmy_sleep(index+1)
                
                start_time_collection = [start_time_collection; filter_start_times(index+1,1)];
                end_time_collection = [end_time_collection; filter_end_times(index+1,1)];
                sleep_collection = [sleep_collection; filter_end_times(index+1,1) - filter_start_times(index+1,1)];
                index = index + 1;
                wcount = wcount + 1;
                
            end
            
            if wcount > 0
                index = index + 1;
            end
            
        end
        
        if length(start_time_collection) < 2
            
            filter_data = [filter_data; {filter_dmy_sleep(j), start_time_collection, end_time_collection, sleep_collection}];
            
        else
            
            wbs = start_time_collection(2:end,:) - end_time_collection(1:length(end_time_collection)-1,:);
            wbs = 24 * wbs;
            %wbs: waking time between sleep episodes (hours)
            
            if min(wbs) >=  wbstime
                
                filter_data = [filter_data; {filter_dmy_sleep(j), start_time_collection, end_time_collection, sleep_collection}];
                
            else
                
                r_start_time_collection = [];
                r_end_time_collection = [];
                r_sleep_collection = [];
                
                frag = wbs < wbstime;
                findex = 1;
                
                for jj = 1:length(frag)
                    
                    if frag(jj) == logical(1) && jj == findex
                        
                        start_time_candidate = start_time_collection(jj);
                        end_time_candidate = end_time_collection(jj+1);
                        
                        if jj < length(frag)
                            
                            fcount = 0;
                            
                            while findex < length(frag) && frag(findex + 1) == logical(1)
                                
                                end_time_candidate = end_time_collection(findex+2);
                                findex = findex + 1;
                                fcount = fcount + 1;
                                
                            end
                            
                            if fcount > 0
                               findex = findex + 1; 
                            end
                            
                        end
                        
                        r_start_time_collection = [r_start_time_collection; start_time_candidate];
                        r_end_time_collection = [r_end_time_collection; end_time_candidate];
                        r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                        
                    else
                        
                        start_time_candidate = start_time_collection(jj);
                        end_time_candidate = end_time_collection(jj);
                        r_start_time_collection = [r_start_time_collection; start_time_candidate];
                        r_end_time_collection = [r_end_time_collection; end_time_candidate] ;
                        r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                        
                        if jj == length(frag)
                            
                            start_time_candidate = start_time_collection(jj+1);
                            end_time_candidate = end_time_collection(jj+1);
                            r_start_time_collection = [r_start_time_collection; start_time_candidate];
                            r_end_time_collection = [r_end_time_collection; end_time_candidate] ;
                            r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                            
                        end
                        
                    end
                    
                    if jj == findex 
                        findex = findex + 1;
                    end
                    
                end
                
                filter_data = [filter_data; {filter_dmy_sleep(j), r_start_time_collection, r_end_time_collection, r_sleep_collection}];
                
            end
            
        end
        
    end
    
    if j == index
        index = index + 1;
    end
    
end

%% Estimation of sleep phase from a sleep file

% 0 - wake, 1 - sleep

Tawake = Tsleep(Tsleep(:,2)==0,:);
dwake = diff(Tawake(:,1));
dsleep_g2 = dwake>(min_sleep_duration/24);

start_times = Tawake(dsleep_g2,1);
end_times = Tawake([logical(0); dsleep_g2],1);

dmy_sleep = floor(start_times);

index = 1;
sleep_data = {};

for j = 1:length(dmy_sleep)
    
    if j == index
        
        start_time_collection = [];
        end_time_collection = [];
        sleep_collection =[];
        
        start_time_collection = [start_time_collection; start_times(j,1)];
        end_time_collection = [end_time_collection; end_times(j,1)];
        sleep_collection = [sleep_collection; end_times(j,1) - start_times(j,1)];
        
        if j < length(dmy_sleep)
            
            wcount = 0;
            
            while index < length(dmy_sleep) && dmy_sleep(j) == dmy_sleep(index+1)
                
                start_time_collection = [start_time_collection; start_times(index+1,1)];
                end_time_collection = [end_time_collection; end_times(index+1,1)];
                sleep_collection = [sleep_collection; end_times(index+1,1) - start_times(index+1,1)];
                index = index + 1;
                wcount = wcount + 1;
                
            end
            
            if wcount > 0
                index = index + 1; 
            end
            
        end
        
        if length(start_time_collection) < 2
            
            sleep_data = [sleep_data; {dmy_sleep(j), start_time_collection, end_time_collection, sleep_collection}];
            
        else
            
            wbs = start_time_collection(2:end,:) - end_time_collection(1:length(end_time_collection)-1,:);
            wbs = 24 * wbs;
            %wbs: waking time between sleep episodes (hours)
            
            if min(wbs) >=  wbstime
                
                sleep_data = [sleep_data; {dmy_sleep(j), start_time_collection, end_time_collection, sleep_collection}];
                
            else
                
                r_start_time_collection = [];
                r_end_time_collection = [];
                r_sleep_collection = [];
                
                frag = wbs < wbstime;
                findex = 1;
                
                for jj = 1:length(frag)
                    
                    if frag(jj) == logical(1) && jj == findex
                        
                        start_time_candidate = start_time_collection(jj);
                        end_time_candidate = end_time_collection(jj+1);
                        
                        if jj < length(frag)
                            
                            fcount = 0;
                            
                            while findex < length(frag) && frag(findex + 1) == logical(1)
                                
                                end_time_candidate = end_time_collection(findex+2);
                                findex = findex + 1;
                                fcount = fcount + 1;
                                
                            end
                            
                            if fcount > 0
                                findex = findex + 1;
                            end
                            
                        end
                        
                        r_start_time_collection = [r_start_time_collection; start_time_candidate];
                        r_end_time_collection = [r_end_time_collection; end_time_candidate];
                        r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                        
                    else
                        
                        start_time_candidate = start_time_collection(jj);
                        end_time_candidate = end_time_collection(jj);
                        r_start_time_collection = [r_start_time_collection; start_time_candidate];
                        r_end_time_collection = [r_end_time_collection; end_time_candidate] ;
                        r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                        
                        if jj == length(frag)
                            
                            start_time_candidate = start_time_collection(jj+1);
                            end_time_candidate = end_time_collection(jj+1);
                            r_start_time_collection = [r_start_time_collection; start_time_candidate];
                            r_end_time_collection = [r_end_time_collection; end_time_candidate] ;
                            r_sleep_collection = [r_sleep_collection; end_time_candidate - start_time_candidate];
                            
                        end
                        
                    end
                    
                    if jj == findex
                       findex = findex + 1; 
                    end
                    
                end
                
                sleep_data = [sleep_data; {dmy_sleep(j), r_start_time_collection, r_end_time_collection, r_sleep_collection}];
                
            end
            
        end
        
    end
    
    if j == index
        index = index + 1;
    end
    
end

%% merge sleep phase data collected from activity file and sleep file

merge_data = [];
sleep_data_size = size(sleep_data);
filter_data_size = size(filter_data);

if isempty(filter_data) == 0
    filter_date = cell2mat(filter_data(:,1));
else
    filter_date = [];
end

filter_index = [];

for k = 1:sleep_data_size(1)
    
    s_date = sleep_data{k,1};
    s_start = sleep_data{k,2};
    s_end = sleep_data{k,3};
    s_sleep = sleep_data{k,4};
    
    if logical(sum(filter_date == s_date)) == logical(0)
        
        [maxsleep,maxindex] = max(s_sleep);
        sleep_onset = mod(s_start(maxindex)*24,24);
        sleep_offset = mod(s_end(maxindex)*24,24);
        mid_point = mod(((s_start(maxindex) + s_end(maxindex))/2)*24,24);
        merge_data =[merge_data; [s_date, mid_point, sleep_onset, sleep_offset]];
        
    else
        
        [maxindex, oindex] = max(filter_date == s_date);
        filter_index = [filter_index; oindex];
        
        f_date = filter_data{oindex,1};
        f_start = filter_data{oindex,2};
        f_end = filter_data{oindex,3};
        f_sleep = filter_data{oindex,4};
        
        m_start = [s_start; f_start];
        m_end =[s_end; f_end];
        m_sleep = [s_sleep; f_sleep];
        m_set = sortrows([m_start, m_end, m_sleep]);
        m_set_size = size(m_set);
        m_difference = m_set(2:end,1) - m_set(1:m_set_size(1)-1,2);
        
        if min(m_difference) > 0
            
            [maxsleep,maxindex] = max(m_set(:,3));
            sleep_onset = mod(m_set(maxindex,1)*24,24);
            sleep_offset = mod(m_set(maxindex,2)*24,24);
            mid_point = mod(((m_set(maxindex,1) + m_set(maxindex,2))/2)*24,24);
            merge_data =[merge_data; [s_date, mid_point, sleep_onset, sleep_offset]];
            
        else
            
            mm_set = [];
            m_index = 1;
            
            for u = 1:length(m_difference)
                
                if m_difference(u) <= 0 && u == m_index
                    
                    m_start_candidate = m_set(u,1);
                    m_end_candidate = max([m_set(u,2), m_set(u+1,2)]);
                    
                    if u < length(m_difference)
                        
                        mcount = 0;
                        
                        while m_index < length(m_difference) && m_difference(m_index+1) <= 0
                            
                            m_end_candidate = max([m_end_candidate, m_set(m_index + 2,2)]);
                            m_index = m_index + 1;
                            mcount = mcount + 1;
                            
                        end
                        
                        if mcount > 0
                            m_index = m_index + 1;
                        end
                        
                    end
                    
                    mm_set = [mm_set; [m_start_candidate, m_end_candidate, m_end_candidate-m_start_candidate]];
                    
                    
                else
                    
                    m_start_candidate =m_set(u,1);
                    m_end_candidate =m_set(u,2);
                    mm_set = [mm_set; [m_start_candidate, m_end_candidate, m_end_candidate-m_start_candidate]];
                    
                end
                
                if u == m_index
                   m_index = m_index +1; 
                end
                               
            end
            
            [maxsleep,maxindex] = max(mm_set(:,3));
            sleep_onset = mod(mm_set(maxindex,1)*24,24);
            sleep_offset = mod(mm_set(maxindex,2)*24,24);
            mid_point = mod(((mm_set(maxindex,1) + mm_set(maxindex,2))/2)*24,24);
            merge_data =[merge_data; [s_date, mid_point, sleep_onset, sleep_offset]];
            
        end
        
    end
    
end

if filter_data_size(1) > 0
    
    filter_index_matrix = 1:filter_data_size(1);
    filter_index_matrix(filter_index) = [];
    
    for kk = 1:length(filter_index_matrix)
        
        f_date = filter_data{filter_index_matrix(kk),1};
        f_start = filter_data{filter_index_matrix(kk),2};
        f_end = filter_data{filter_index_matrix(kk),3};
        f_sleep = filter_data{filter_index_matrix(kk),4};
        
        [maxsleep,maxindex] = max(f_sleep);
        sleep_onset = mod(f_start(maxindex)*24,24);
        sleep_offset = mod(f_end(maxindex)*24,24);
        mid_point = mod(((f_start(maxindex) + f_end(maxindex))/2)*24,24);
        
        merge_data =[merge_data; [f_date, mid_point, sleep_onset, sleep_offset]];
        
    end
    
end

merge_data = sortrows(merge_data);
merge_size = size(merge_data);

for kkk =1:merge_size(1)
    
    fprintf(results_file, '%s, %f, %f, %f, %f\n', datestr(merge_data(kkk,1)), merge_data(kkk,2), merge_data(kkk,3), merge_data(kkk,4), exist_time_zone);
    
end

end


