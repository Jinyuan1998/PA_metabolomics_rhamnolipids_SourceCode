close all; clear;

addpath('../matlabfunction');

tsdpath = 'tsd_glycerol';

fileformat = fullfile(tsdpath, '*.mat')
tsdfiles = dir(fileformat);
tsdfiles = {tsdfiles.name};

for i=1:length(tsdfiles)

    t = load([tsdpath '/' tsdfiles{i}]);
    t = t.tsd;
    if i==1
        T=t;
    else
        T = [T t];
    end
    clear t
end
%%
% cmap = [.6 .6 .9; 
%         .2 .8 .6;
%         .8 .5 .5;
%         0.4    0   0.6;
%         .2  0 .5];
lw = 1;
lwd = 0;
% cmap = [
%     0  255 255;
%     1 1 1;
%      
%      50 154 250;
%      220 61 160;
%      50 190 120;
%      115 82 68;
%      195 150 130;
%      8 133 170 ] / 255;
% cmap = hsv(5)*0.9;
% cmap(end, :) =[150 140 151] / 255;
% % cmap(end, :) =[215 40 251] / 255;
% cmap(2, :) =[1 1  1] / 255;
cmap = [
     255 103 70;
    2,2,2;  
    255 215 0;    
     0,191,255 ;
     155,100,154;
     ] / 255;
alpha=0;    
od_range = [0.1 1.6];
gfp_range = [ 100 2e4];
strains = {'PA14-con' 'PA14-PrhlAB'};
stress = {'H2O2-0' 'H2O2-10mM' 'H2O2-20mM' 'H2O2-50mM' 'H2O2-100mM'};
indeces = repmat('[]', length(strains), length(stress));
figure
subplot(2,2,1)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));

    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             od_norm = (T(idx(k)).wavelength(1).data) .^ 1.1027* exp(0.0069);
            plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        else
        plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        hold on
    end
end
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, PA1/04/03', 'fontsize', 14)

subplot(2,2,2)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
            plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
        end
        hold on
    end
end
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, PA1/04/03', 'fontsize', 14)


subplot(2,2,3)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
            plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
        plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        hold on
    end
end
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, PrhlAB', 'fontsize', 14)

subplot(2,2,4)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             gfp_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(2).data,  normTbl.time);
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
            plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*j], 'linewidth', lw-lwd*j);        
        else
            pl=plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
%             pl.Color(4) = 0.2;
        end
        hold on
    end
end
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, PrhlAB', 'fontsize', 14)
set(gcf, 'Renderer', 'painters')
%%
figure
for i=1:length(stress)
    plot([1 2], [i i], '-', 'Color', cmap(i,:),'linewidth', 5)
    hold on
    
end
legend(stress, 'Location', 'southoutside')

%% plot mean and range
lw1=1;
alp = 0.3;
time_delta = 0.5;
figure
subplot(2,2,1)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));    
    for k=1:length(idx)
        addRO1 = T(idx(k)).addRO1;
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             od_norm = (T(idx(k)).wavelength(1).data) .^ 1.1027* exp(0.0069);

%             plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        else
            od_norm = T(idx(k)).wavelength(1).data;
            
        end
        od_norm_mean = mean(od_norm,2);
        od_norm_max = max(od_norm, [],2);
        od_norm_min = min(od_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                 [od_norm_min(Idx), od_norm_max(Idx)-od_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), od_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on

line([addRO1 addRO1], od_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, PA14::PA1/04/03-GFP', 'fontsize', 14)

subplot(2,2,2)
for j=1:length(stress)
    
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));
    for k=1:length(idx)
        addRO1 = T(idx(k)).addRO1;
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
%             plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            gfp_norm = T(idx(k)).wavelength(2).data;
%             plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
        end
        gfp_norm_mean = mean(gfp_norm,2);
        gfp_norm_max = max(gfp_norm, [],2);
        gfp_norm_min = min(gfp_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                [gfp_norm_min(Idx), gfp_norm_max(Idx)-gfp_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        plot(T(idx(k)).time(Idx), gfp_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on

    end
end
hold on
line([addRO1 addRO1], gfp_range, 'LineWidth', 2, 'color', 'k')        
hold off
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, PA14::PA1/04/03-GFP', 'fontsize', 14)


subplot(2,2,3)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));    
    for k=1:length(idx)
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            od_norm = T(idx(k)).wavelength(1).data;
%         plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        od_norm_mean = mean(od_norm,2);
        od_norm_max = max(od_norm, [],2);
        od_norm_min = min(od_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                [od_norm_min(Idx), od_norm_max(Idx)-od_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), od_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on
% addRO1 = T(idx(k)).addRO1;
line([addRO1 addRO1], od_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, PA14::PrhlAB-GFP', 'fontsize', 14)

subplot(2,2,4)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        Idx = find(T(idx(k)).time >= addRO1 + time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             gfp_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(2).data,  normTbl.time);
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
%             plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*j], 'linewidth', lw-lwd*j);        
        else
            gfp_norm = T(idx(k)).wavelength(2).data;
%             pl=plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
%             pl.Color(4) = 0.2;
        end
        gfp_norm_mean = mean(gfp_norm,2);
        gfp_norm_max = max(gfp_norm, [],2);
        gfp_norm_min = min(gfp_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), [gfp_norm_min(Idx), gfp_norm_max(Idx)-gfp_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), gfp_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on
% % addRO1 = T(idx(k)).addRO1;
line([addRO1 addRO1], gfp_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, PA14::PrhlAB-GFP', 'fontsize', 14)
set(gcf, 'Renderer', 'painters')
%%
%%

% cmap = [.6 .6 .9; 
%         .2 .8 .6;
%         .8 .5 .5;
%         0.4    0   0.6;
%         .2  0 .5];
lw = 1;
lwd = 0;
% cmap = [
%     0  255 255;
%     1 1 1;
%      
%      50 154 250;
%      220 61 160;
%      50 190 120;
%      115 82 68;
%      195 150 130;
%      8 133 170 ] / 255;
% cmap = hsv(5)*0.9;
% cmap(end, :) =[150 140 151] / 255;
% % cmap(end, :) =[215 40 251] / 255;
% cmap(2, :) =[1 1  1] / 255;
cmap =[
     255 103 70;
    2,2,2;  
    255 215 0;    
     0,191,255 ;
     155,100,154;
     ] / 255 ;
alpha=0;    
od_range = [0.1 1.6];
gfp_range = [ 100 3e4];
strains = {'rhlA-con' 'rhlA-PrhlAB'};
stress = {'H2O2-0' 'H2O2-10mM' 'H2O2-20mM' 'H2O2-50mM' 'H2O2-100mM'};
indeces = repmat('[]', length(strains), length(stress));
%%
figure
subplot(2,2,1)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));

    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             od_norm = (T(idx(k)).wavelength(1).data) .^ 1.1027* exp(0.0069);
            plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        else
        plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        hold on
    end
end
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, -rhlA::PA1/04/03-GFP', 'fontsize', 14)

subplot(2,2,2)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
            plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
        end
        hold on
    end
end
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, -rhlA::PA1/04/03-GFP', 'fontsize', 14)


subplot(2,2,3)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
            plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
        plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        hold on
    end
end
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, -rhlA::PrhlAB-GFP', 'fontsize', 14)

subplot(2,2,4)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        if strcmp(T(idx(k)).tecan, 'spark')
%             gfp_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(2).data,  normTbl.time);
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
            plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*j], 'linewidth', lw-lwd*j);        
        else
            pl=plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
%             pl.Color(4) = 0.2;
        end
        hold on
    end
end
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, -rhlA::PrhlAB-GFP', 'fontsize', 14)
set(gcf, 'Renderer', 'painters')
%% plot mean and range
lw1=2;
alp = 0.3;
time_delta = 0.5;
figure
subplot(2,2,1)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));    
    for k=1:length(idx)
        addRO1 = T(idx(k)).addRO1;
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             od_norm = (T(idx(k)).wavelength(1).data) .^ 1.1027* exp(0.0069);

%             plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        else
            od_norm = T(idx(k)).wavelength(1).data;
            
        end
        od_norm_mean = mean(od_norm,2);
        od_norm_max = max(od_norm, [],2);
        od_norm_min = min(od_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                 [od_norm_min(Idx), od_norm_max(Idx)-od_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), od_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on

line([addRO1 addRO1], od_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, -rhlA::PA1/04/03-GFP', 'fontsize', 14)

subplot(2,2,2)
for j=1:length(stress)
    
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{1}));
    for k=1:length(idx)
        addRO1 = T(idx(k)).addRO1;
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
%             plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            gfp_norm = T(idx(k)).wavelength(2).data;
%             plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
        end
        gfp_norm_mean = mean(gfp_norm,2);
        gfp_norm_max = max(gfp_norm, [],2);
        gfp_norm_min = min(gfp_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                [gfp_norm_min(Idx), gfp_norm_max(Idx)-gfp_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        plot(T(idx(k)).time(Idx), gfp_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on

    end
end
hold on
line([addRO1 addRO1], gfp_range, 'LineWidth', 2, 'color', 'k')        
hold off
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, -rhlA::PA1/04/03-GFP', 'fontsize', 14)


subplot(2,2,3)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));    
    for k=1:length(idx)
        Idx = find(T(idx(k)).time >= addRO1 +time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             od_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(1).data,  normTbl.time);
            od_norm = T(idx(k)).wavelength(1).data * 0.7291 + 0.0033;
%             plot(T(idx(k)).time, od_norm, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j);
        
        else
            od_norm = T(idx(k)).wavelength(1).data;
%         plot(T(idx(k)).time, T(idx(k)).wavelength(1).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); %, 'MarkerEdgeColor', cmap(j,:))
        end
        od_norm_mean = mean(od_norm,2);
        od_norm_max = max(od_norm, [],2);
        od_norm_min = min(od_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), ...
                [od_norm_min(Idx), od_norm_max(Idx)-od_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), od_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on
% addRO1 = T(idx(k)).addRO1;
line([addRO1 addRO1], od_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', od_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('OD600 [AU]', 'fontsize', 12)
title('OD600, -rhlA::PrhlAB-GFP', 'fontsize', 14)

subplot(2,2,4)
for j=1:length(stress)
    idx = find(contains({T.name}, stress(j)) & contains({T.name}, strains{2}));
    for k=1:length(idx)
        Idx = find(T(idx(k)).time >= addRO1 + time_delta);
        if strcmp(T(idx(k)).tecan, 'spark')
%             gfp_intpl = interp1(T(idx(k)).time, T(idx(k)).wavelength(2).data,  normTbl.time);
            gfp_norm = T(idx(k)).wavelength(2).data * 0.5119-19.483;
%             plot(T(idx(k)).time, gfp_norm, '-', 'Color', [cmap(j,:) 1-alpha*j], 'linewidth', lw-lwd*j);        
        else
            gfp_norm = T(idx(k)).wavelength(2).data;
%             pl=plot(T(idx(k)).time, T(idx(k)).wavelength(2).data, '-', 'Color', [cmap(j,:) 1-alpha*(j-1)], 'linewidth', lw-lwd*j); 
%             pl.Color(4) = 0.2;
        end
        gfp_norm_mean = mean(gfp_norm,2);
        gfp_norm_max = max(gfp_norm, [],2);
        gfp_norm_min = min(gfp_norm, [], 2);        
        h = area(T(idx(k)).time(Idx), [gfp_norm_min(Idx), gfp_norm_max(Idx)-gfp_norm_min(Idx)], 'FaceColor', cmap(j,:));
        set(h, 'EdgeColor', 'none');
        % erase the first area
        set(h(1), 'FaceColor', 'none');
        h(2).FaceAlpha = alp;
        hold on
        plot(T(idx(k)).time(Idx), gfp_norm_mean(Idx), '-', 'Color', [cmap(j,:) 1], 'linewidth', lw1); %, 'MarkerEdgeColor', cmap(j,:))
        hold on
    end
end
hold on
% % addRO1 = T(idx(k)).addRO1;
line([addRO1 addRO1], gfp_range, 'LineWidth', 2, 'color', 'k')
hold off
set(gca, 'Ylim', gfp_range)
set(gca, 'YScale', 'log')
xlabel('time (h)', 'fontsize', 12)
ylabel('GFP [AU]', 'fontsize', 12)
title('GFP, -rhlA::PrhlAB-GFP', 'fontsize', 14)
set(gcf, 'Renderer', 'painters')