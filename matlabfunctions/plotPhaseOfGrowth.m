%% plot growth phases
function plotPhaseOfGrowth(grwothCurve, timepoints, phase_start_index, color_backbone)
% figure();
for i=1:size(grwothCurve,2)
%     subplot(5,7,i);
%     hold on;
    if nargin==3
        bColor = [.5 .5 .5];
    else
        bColor =  color_backbone;
    end
    yobs = grwothCurve(:,i);  
    
    lw2 = 5;
    dgealpha2=.3;
    if (sum(isnan(phase_start_index(i,:)))==0)
        patchline(timepoints(phase_start_index(i,1):phase_start_index(i,2)), yobs(phase_start_index(i,1):phase_start_index(i,2)),...
            'linestyle', '-', 'edgecolor', [1 1 0], 'linewidth', lw2, 'edgealpha', 0.5);
        patchline(timepoints(phase_start_index(i,2):phase_start_index(i,3)), yobs(phase_start_index(i,2):phase_start_index(i,3)),...
            'linestyle', '-', 'edgecolor', [.5 1 1], 'linewidth', lw2, 'edgealpha', 0.5);
%         patchline(timepoints(phase_start_index(i,3):end), yobs(phase_start_index(i,3):end),...
%             'linestyle', '-', 'edgecolor', [1 .2 .8], 'linewidth', lw2, 'edgealpha', 0.2);
        patchline(timepoints(phase_start_index(i,3):end), yobs(phase_start_index(i,3):end),...
            'linestyle', '-', 'edgecolor', [.8 .8 .8], 'linewidth', lw2, 'edgealpha', 0.5);
    end
    hold on
    plot(timepoints, yobs, '-', 'LineWidth', 2, 'color', bColor);
    
%     if (sum(isnan(phase_start_index(i,:)))==0)
%         patchline(timepoints(phase_start_index(i,1):phase_start_index(i,2)), yobs(phase_start_index(i,1):phase_start_index(i,2)),...
%             'linestyle', '-', 'edgecolor', [.1 .1 .1], 'linewidth', lw2, 'edgealpha', dgealpha2);
%         patchline(timepoints(phase_start_index(i,2):phase_start_index(i,3)), yobs(phase_start_index(i,2):phase_start_index(i,3)),...
%             'linestyle', '-', 'edgecolor', [.6 .6 .6], 'linewidth', lw2, 'edgealpha', dgealpha2);
% %         patchline(timepoints(phase_start_index(i,3):phase_start_index(i,end)), yobs(phase_start_index(i,3):phase_start_index(i,end)),...
% %             'linestyle', '-', 'edgecolor', [.9 .9 .9], 'linewidth', lw2, 'edgealpha', dgealpha2);
%     end


%     xlim([0,48]);
%     set(gca,'XTick',[0,12,24,36,48]);
%     ylim([0,1.5]);
%     set(gca,'YTick',[0,0.5,1.0,1.5]);
%     axis square;
    box on;
    xlabel('Time (hour)');
%     ylabel('OD');
%     title(titles);
end

%% save to file
% tbl_phase_start_time = array2table(phase_start_time);
% tbl_phase_start_time.Properties.RowNames = tblgc.Properties.VariableNames;
% tbl_phase_start_time.Properties.VariableNames = {'Phase1';'Phase2';'Phase3'};
% writetable(tbl_phase_start_time, 'PA_glycerol_growth_phase_start_time.csv', 'Delimiter', ',', 'WriteRowNames', true);