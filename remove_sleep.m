function new_Thr = remove_sleep(Thr, Tsleep)

%%% Function that takes a heart rate and sleep array and removes any heart
%%% rate data that overlaps with sleep data 

%% Account for the first part of the heart rate data

v = find(Tsleep(:,2)>0);

new_Thr = Thr(Thr(:,1)<Tsleep(v(1),1),:);

%% Take the diff

Tasleep = Tsleep(Tsleep(:,2)>0,:);

dsleep = diff(Tasleep(:,1));

dsleep_g2 = dsleep>(2/24);

start_times = Tasleep([logical(0); dsleep_g2],1);
end_times = Tasleep(dsleep_g2,1);

wake_periods = [end_times start_times];

for i = 1:length(wake_periods(:,1))

    new_Thr = [new_Thr; Thr(Thr(:,1)>wake_periods(i,1) & Thr(:,1)<wake_periods(i,2),:)];

end

new_Thr = [new_Thr; Thr(Thr(:,1)>Tsleep(v(end),1),:)];

%figure(1)
%plot(Thr(:,1), Thr(:,2))
%hold on
%plot(new_Thr(:,1), new_Thr(:,2))
%hold off


