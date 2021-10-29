function growthRate=extractGrowthRate(d1, times, phase_start_time)

growthRate = zeros(size(d1, 2), 2); % first column is median growth rate in phase I, second is phase II
for j=1:size(d1, 2)
    od_delta = ppval(fnder(spline(times, d1(:,j)),1), times); % first order derivative
    od_delta_smoothed = sgolayfilt(od_delta, 3, 51); % smoothing first-order derivative
    d1_smooth=sgolayfilt(d1(:,j), 3, 51);
    growthrate = od_delta_smoothed' ./d1_smooth;
    growthRate(j, 1) = median(growthrate(times>=phase_start_time(j,1) & times<= phase_start_time(j,2)));
    growthRate(j, 2) = median(growthrate(times>phase_start_time(j,2) & times<= phase_start_time(j,3)));
end