function sleep_filter(Tsleep, min_sleep_duration, wbstime, exist_time_zone)

%%% Function that calculates the midpoint of sleep.

% Open results.csv file and add a header line
results_file = fopen('Sleep_results.csv', 'w');
fprintf(results_file, 'Date (dd-mmm-yyyy), sleep midpoint, sleep onset, sleep offset, Existence of Time Zone Information \n');

% 0 - wake, 1 - sleep

Tawake = Tsleep(Tsleep(:,2)==0,:);
dwake = diff(Tawake(:,1));
dsleep_g2 = dwake>(min_sleep_duration/24);

start_times = Tawake(dsleep_g2,1);
end_times = Tawake([logical(0); dsleep_g2],1);

%[sleep_days, ~, idx] = unique(floor(start_times), 'stable');    
%dmy_sleep = sleep_days;

dmy_sleep = floor(start_times);

index = 1;

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
            
            sleep_onset = mod(start_time_collection(1)*24,24);
            sleep_offset = mod(end_time_collection(1)*24,24);
            mid_point = mod(((start_time_collection(1) + end_time_collection(1))/2)*24,24);
            fprintf(results_file, '%s, %f, %f, %f, %f\n', datestr(dmy_sleep(j)), mid_point, sleep_onset, sleep_offset, exist_time_zone);
            
        else
            
            wbs = start_time_collection(2:end,:) - end_time_collection(1:length(end_time_collection)-1,:);
            wbs = 24 * wbs;
            %wbs: waking time between sleep episodes (hours)
            
            if min(wbs) >=  wbstime
                
                [maxsleep, maxindex] = max(sleep_collection);
                sleep_onset = mod(start_time_collection(maxindex)*24,24);
                sleep_offset = mod(end_time_collection(maxindex)*24,24);
                mid_point = mod(((start_time_collection(maxindex) + end_time_collection(maxindex))/2)*24,24);
                fprintf(results_file, '%s, %f, %f, %f, %f\n', datestr(dmy_sleep(j)), mid_point, sleep_onset, sleep_offset, exist_time_zone);
                
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
                
                [maxsleep, maxindex] = max(r_sleep_collection);
                sleep_onset = mod(r_start_time_collection(maxindex)*24,24);
                sleep_offset = mod(r_end_time_collection(maxindex)*24,24);
                mid_point = mod(((r_start_time_collection(maxindex) + r_end_time_collection(maxindex))/2)*24,24);
                fprintf(results_file, '%s, %f, %f, %f, %f\n', datestr(dmy_sleep(j)), mid_point, sleep_onset, sleep_offset, exist_time_zone);
                
            end
            
        end
        
    end
    
    if j == index
        index = index + 1; 
    end

end


end