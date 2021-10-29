% This script divides each growth curve into phase I, II and III
% last modified by Chen Liao, May 03, 2021
function [phase_start_index, phase_start_time] = phaseExtraction(growthCurve, timepoints, varargin)
%% read growth curves
% tblgc = readtable('normalized_mean_growth_curve_PA_glycerol.csv');
% tblgc=growthCurveT;
% tblgc.Properties.VariableNames{1} = 'Time';
% timepoints = tblgc{:,1};
% ntps = length(timepoints);
% tblgc.Time = [];
if isempty(varargin)
    needSmooth = 0;
    carbon='glycerol';
elseif isempty(varargin{1})
    needSmooth = 0;
elseif length(varargin) == 2 & isempty(varargin{2}) 
    carbon='glycerol';
elseif nargin > 4
    disp('too many inputs!!')
elseif nargin == 4
    needSmooth = varargin{1};
    carbon=varargin{2};
end
if strcmp(carbon, 'glycerol')    
    yobs_prime_smooth_cutoff=0.001;
elseif strcmp(carbon, 'glucose')
    yobs_prime_smooth_cutoff=0.02;
end
%% identify growth phases
phase_start_index = zeros(size(growthCurve,2),3); % number of growth curves x 3
phase_start_time = zeros(size(growthCurve,2),3); % number of growth curves x 3

% loop over each growth curve
for i=1:size(growthCurve,2)
    
    % calculate first-order derivative
    if needSmooth == 1 
        yobs = sgolayfilt(growthCurve(:, i),3, 51); % smooth observed data
    else
        yobs = growthCurve(:, i);
    end
    yobs_prime = ppval(fnder(spline(timepoints, yobs),1), timepoints); % first order derivative
    yobs_prime_smoothed = sgolayfilt(yobs_prime, 3, 51); % smoothing first-order derivative
    
    % find growth phases
%     index_pos_1st_deriv = find(yobs_prime_smoothed>0.001); % index of timepoints where first derivatives are larger than a threshold
    index_pos_1st_deriv = find(yobs_prime_smoothed> yobs_prime_smooth_cutoff);    
    if (~isempty(index_pos_1st_deriv))       
        % Phase I & II: find the longest interval that first derivative is above 0
        %
        % It is possible that there are temporary fluctuations where first
        % derivative drops below 0 for a short time period and comes back.
        % Ignore these temporary drops using a parameter "max_gap_length"
        max_gap_length = 10;
        
        % maximum length of phase 12 found throughout the entire region of
        % positive derivatives
        max_interval_length = 0;
        
        % the start index of phase 12 and phase 3 for the current region of
        % positive derivatives (gaps less than max_gap_length are ignored)
        curr_phase_12_start_index = index_pos_1st_deriv(1);
        curr_phase_3_start_index = index_pos_1st_deriv(1);
        
        % the start index of phase 12 and phase 3 after searching for the
        % entire region of positive derivatives
        final_phase12_start_index = curr_phase_12_start_index;
        final_phase3_start_index = curr_phase_3_start_index;
        
        for j=2:length(index_pos_1st_deriv)
            if ((index_pos_1st_deriv(j)-index_pos_1st_deriv(j-1))<=max_gap_length)
                % as long as the gap is less than max_gap_length, ignore
                % the gap
                curr_phase_3_start_index = index_pos_1st_deriv(j); % move the right end to current index
                % if the last index and the current phase is longer than
                % any phase found ever, then update the phase to the
                % current phase
                if (j==length(index_pos_1st_deriv) & curr_phase_3_start_index-curr_phase_12_start_index > max_interval_length)
                    final_phase12_start_index = curr_phase_12_start_index;
                    final_phase3_start_index = curr_phase_3_start_index;
                end
            else
                % if the gap is longer than max_gap_length, then update the
                % phase information (only update when the current phase is
                % longer than max_interval_length) and move to the next
                % region
                if (curr_phase_3_start_index-curr_phase_12_start_index > max_interval_length)
                    final_phase12_start_index = curr_phase_12_start_index;
                    final_phase3_start_index = curr_phase_3_start_index;
                    max_interval_length = final_phase3_start_index - final_phase12_start_index;
                end
                curr_phase_12_start_index = index_pos_1st_deriv(j);
                curr_phase_3_start_index = index_pos_1st_deriv(j);
            end
        end
        phase_start_index(i,1) = final_phase12_start_index;
        phase_start_index(i,3) = final_phase3_start_index;
        phase_start_time(i,1) = timepoints(final_phase12_start_index);
        phase_start_time(i,3) = timepoints(final_phase3_start_index);
        
        % Find separation of phase 1 and phase 2
        % Phase II: if multiple peaks exist within the interval, choose the
        % highest peak as the entry point for Phase II
        %
        % find all peaks within the region of positive first order
        % derivative
%         [x_max, y_max, A]=crit_interp_g(yobs_prime_smoothed, timepoints);
%         xq=timepoints(1):0.01:timepoints(end);
%         V = interp1( timepoints, yobs_prime_smoothed,xq, 'spline');
%         figure; plot(xq, V, 'bo', timepoints, yobs_prime_smoothed, 'gx', timepoints, growthCurve/5, 'k.', timepoints(locs), pks, 'm*')
%         
        [pks, locs] = findpeaks(yobs_prime_smoothed);
       
%         figure; plot(timepoints, yobs_prime_smoothed, 'go', timepoints(locs), yobs_prime_smoothed(locs),'rx')
        % filter out the peaks located outside the phase 12 region
        pks2 = pks(find(locs >= final_phase12_start_index & locs <= final_phase3_start_index));
        locs2 = locs(find(locs >= final_phase12_start_index & locs < final_phase3_start_index));
        if (~isempty(locs2))
            [~,max_peak_index] = max(pks2);
            phase_start_index(i,2) = locs2(max_peak_index);
            phase_start_time(i,2) = timepoints(locs2(max_peak_index));
        elseif isempty(locs2) && ~isempty(locs)
%             [~,max_peak_index] = max(pks);
%             idx = locs(max_peak_index);
            a = diff(yobs_prime_smoothed(1:(final_phase3_start_index-1)));
            
            if a <0
                phase_start_index(i,2)=phase_start_index(i,1)+1;
                phase_start_time(i,2) = timepoints(phase_start_index(i,2));
            end
        else
            phase_start_index(i,2) = NaN;
            phase_start_time(i,2) = NaN;
        end
    else
        % first-order derivative is smaller than a threshold for the entire
        % time course.
        fprintf('%s: no growth at all.\n', tblgc.Properties.VariableNames{i});
        phase_start_index(i,:) = [NaN, NaN, NaN];
    end
end

% %% plot growth phases
% figure();
% for i=1:size(tblgc,2)
%     subplot(5,7,i);
%     hold on;
%     
%     yobs = tblgc{:,i};  
%     plot(timepoints, yobs, 'k-', 'LineWidth', 1);
%     
%     if (sum(isnan(phase_start_index(i,:)))==0)
%         patchline(timepoints(phase_start_index(i,1):phase_start_index(i,2)), yobs(phase_start_index(i,1):phase_start_index(i,2)),...
%             'linestyle', '-', 'edgecolor', 'r', 'linewidth', 4, 'edgealpha', 0.2);
%         patchline(timepoints(phase_start_index(i,2):phase_start_index(i,3)), yobs(phase_start_index(i,2):phase_start_index(i,3)),...
%             'linestyle', '-', 'edgecolor', 'b', 'linewidth', 4, 'edgealpha', 0.2);
%         patchline(timepoints(phase_start_index(i,3):end), yobs(phase_start_index(i,3):end),...
%             'linestyle', '-', 'edgecolor', 'g', 'linewidth', 4, 'edgealpha', 0.2);
%     end
%     
%     xlim([0,48]);
%     set(gca,'XTick',[0,12,24,36,48]);
%     ylim([0,1.5]);
%     set(gca,'YTick',[0,0.5,1.0,1.5]);
%     axis square;
%     box on;
%     xlabel('Time (hour)');
%     ylabel('OD');
%     title(tblgc.Properties.VariableNames{i});
% end
% 
% %% save to file
% tbl_phase_start_time = array2table(phase_start_time);
% tbl_phase_start_time.Properties.RowNames = tblgc.Properties.VariableNames;
% tbl_phase_start_time.Properties.VariableNames = {'Phase1';'Phase2';'Phase3'};
% writetable(tbl_phase_start_time, 'PA_glycerol_growth_phase_start_time.csv', 'Delimiter', ',', 'WriteRowNames', true);