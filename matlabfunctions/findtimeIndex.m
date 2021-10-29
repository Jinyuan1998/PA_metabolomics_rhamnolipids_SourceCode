function timeIdex = findtimeIndex(timepoints, phase_time)

% phase_time = phase_start_time_median(:,1);
timeIdex = zeros(length(phase_time), 1);
for i=1:length(phase_time)
timeIdex(i)= min(find(timepoints>=phase_time(i)));
end