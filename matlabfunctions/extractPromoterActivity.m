function promoterActivity=extractPromoterActivity(d1, d2, times, phase_start_time) %d1 is the od and d2 is gfp

promoterActivity = zeros(size(d1, 2), 6); % first column is median growth rate in phase I, second is phase II
for j=1:size(d1, 2)
    gfp_delta = ppval(fnder(spline(times, d2(:,j)),1), times); % first order derivative
    gfp_delta_smoothed = sgolayfilt(gfp_delta, 3, 51); % smoothing first-order derivative
    d1_smooth=sgolayfilt(d1(:,j), 3, 51);
    pd = gfp_delta_smoothed' ./ d1_smooth;
    promoterActivity(j, 1) = median(pd(times>=phase_start_time(j,1) & times<= phase_start_time(j,2)));
    promoterActivity(j, 2) = median(pd(times>phase_start_time(j,2) & times<= phase_start_time(j,3)));
    promoterActivity(j, 3) = prctile(pd(times>=phase_start_time(j,1) & times<= phase_start_time(j,2)), 25);
    promoterActivity(j, 4) = prctile(pd(times>phase_start_time(j,2) & times<= phase_start_time(j,3)), 25);
    promoterActivity(j, 5) = prctile(pd(times>=phase_start_time(j,1) & times<= phase_start_time(j,2)), 75);
    promoterActivity(j, 6) = prctile(pd(times>phase_start_time(j,2) & times<= phase_start_time(j,3)), 75);
end

promoterActivity=array2table(promoterActivity, 'VariableNames', ...
                        {'phaseI_median' 'phaseII_median' 'phaseI_25pctl' 'phaseII_25pctl' 'phaseI_75pctl' 'phaseII_75pctl'});