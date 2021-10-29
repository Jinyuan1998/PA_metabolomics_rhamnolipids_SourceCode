classdef TimeSeriesData_jy
    % defines objects that read a Tecan file and interprete the
    % time series data
    
    % properties
    properties (SetAccess = public)
        numberOfWavelengths = 0;
        numberOfTimePoints = 0;
        wavelengthData = [];
        samples = [];
        map = [];
    end
    
    % Method definisions
    methods
%         % constructor: creates a new instance of time series
%         % reads all data from the matrices
%         function tsd = TimeSeriesData(filename)
%             fid = fopen(filename);
%             
%             tline = fgetl(fid);
%             while ischar(tline)
%                 % get the number of time points
%                 firstCell = TimeSeriesData.textInCell(tline, 1);
%                 if ~isempty(strfind(firstCell, 'Cycles'))
%                     tsd.numberOfTimePoints =...
%                         num2str(TimeSeriesData.textInCell(tline, 5));
%                 end;
%                 
%                 % get the data matrices
%                 if ~isempty(strfind(firstCell, 'Cycle Nr.'))
%                     % increment the wavelength count
%                     tsd.numberOfWavelengths = tsd.numberOfWavelengths + 1;
%                     % get the wavelength name from previous line
%                     tsd.wavelengthData(tsd.numberOfWavelengths).name =...
%                         TimeSeriesData.textInCell(previousLine, 1);
%                     % extract well map
%                     commas = strfind(tline, ','); % find indexes of commas
%                     map = [];
%                     column = 1;
%                     for j = 4:length(commas)
%                         map{end+1} = tline(commas(j-1)+1:commas(j)-1);
%                         column = column + 1;
%                     end;
%                     % get the cells in last columns
%                     map{column} = tline(commas(j)+1:end);
%                     tsd.map = map;
%                     % jump to next line
%                     tline = fgetl(fid);
%                     % start reading and compiling data until finding
%                     % empty cell
%                     n = 0;
%                     while TimeSeriesData.numberInCell(tline, 1) == n+1
%                         %check if this line has any data
%                         if TimeSeriesData.checkIfLineHasData(tline),
%                             % extract time, convert to hour
%                             tsd.wavelengthData(tsd.numberOfWavelengths)...
%                                 .times(n+1) =...
%                                 TimeSeriesData.numberInCell(tline, 2) / 3600;
%                             % extract data
%                             try 
%                                 tsd.wavelengthData(tsd.numberOfWavelengths)...
%                                     .data(n+1, :) =...
%                                     TimeSeriesData.numberInLine(tline, 4);
%                             catch
%                                 lixo = 0;
%                             end
%                             % update line number
%                             n = TimeSeriesData.numberInCell(tline, 1);
%                         end;
%                         % proceed to next line in spreadsheet
%                         tline = fgetl(fid);
%                     end;
%                 end;
%                 % proceed to next line
%                 previousLine = tline;
%                 tline = fgetl(fid);
%             end
%             
%             fclose(fid);
%         end
        

% this function read from xlsx file instead of csv file. Saves a step to
% save the xlsx file to csv file
            function tsd = TimeSeriesData_jy(filename, model)
            if nargin == 1
                model = 'other';
            end
            [~, ~, All] = xlsread(filename);
            
            for i=1:size(All, 1)
                B = All(i, :);
                
                % get number of wavelength
%                 x = strfind(B{1},'Labels');
%                 if ~isempty(x)
% %                     noftimepoints =B{5};
%                     wn = B{1};
%                     tsd.numberOfWavelengths = wn(1);                    
%                 end
                
                % get the number of time points
                x = strfind(lower(B{1}), 'cycles');
                if ~isempty(x)
%                     noftimepoints =B{5};
%                     tsd.numberOfTimePoints = B{5};   
                     p = All(i:end, 1);
                     
                     pp = cellfun(@isnan,p,'UniformOutput', false);
                     pp = cellfun(@(X) X(1), pp, 'UniformOutput', false);
                     pp = [pp{:}];
                     pp2 = p(pp ~= 1);
                     if strcmp(model, 'spark')
                        ppp = cellfun(@(X) regexpi(X,'^[0-9]'), pp2 , 'UniformOutput', false);
                        ppp = cellfun(@isempty , ppp );
                        pp3 = pp2(ppp ==0);
                        tsd.numberOfTimePoints = max(cellfun(@str2num, pp3));
                     else
                        ppp = cellfun(@(X) isa(X, 'double'), pp2 );
                        pp3 = pp2(ppp == 1);
                        tsd.numberOfTimePoints = max(cellfun(@(X) X(1), pp3));
                     end
                end
                % get the data matrices
                y = strfind(lower(B{1}), 'cycle nr.');
                if ~isempty(y)
                    tsd.numberOfWavelengths = tsd.numberOfWavelengths + 1;
                    tsd.wavelengthData(tsd.numberOfWavelengths).name = All{i-1,1};
                    if ~isempty(strfind(B{2}, 'Time [s]'))
                        time = All(i+1:i+tsd.numberOfTimePoints, 2);
%                         time = All(i+1:i+2, 2);
                        time = arrayfun(@(X) X{1}, time);
                        tsd.wavelengthData(tsd.numberOfWavelengths).times = time' / 3600;
                    end
                    % get the wells
                    map = All(i,4:end);
                    tsd.map = map;
                    % get the column number
                    column = length(map);
                    % get the data for each wavelength
                    x = All(i+1:i+tsd.numberOfTimePoints, 4:4+column-1);
                    
                    %%% deal with data that could be 'over' the detection
                    for u=1:size(x,2)
                        uu = x(:, u);
                        overidx = find(strcmp(uu, 'OVER'));
                        for uui = 1:length(overidx)
                            uu(overidx(uui)) = uu(overidx(uui)-1); % use the previous data value to fill in the over reached points
                        end
                        x(:, u) = uu;
                    end
                    %%% deal with NaN data
                    for u=1:size(x, 2)
                        uu = x(:, u);
%                         idx = find(isnan(uu));
                        idx = find(strcmp(uu, 'NaN'));
                        for uui = 1:length(idx)
                            uu(overidx(uui)) = uu(overidx(uui)-1);
                        end
                    end
                    
                    
                    xsize = size(x);
                    xx = [x{:}];
                    newx = reshape(xx, xsize);
                    tsd.wavelengthData(tsd.numberOfWavelengths).data = newx;
                    clear x xsize xx newx
                    i = i + tsd.numberOfTimePoints;
                end            
            end
        end

    % Adds a second spreadsheet to this TimeSeriesData.
        % This can be used to join to consecutive time series.
        % Must be added before the samples are set.
        function tsd = appendMoreTimeSeriesData(tsd, filename, gap, model)
            if nargin == 2
                gap = 0;
                model = 'other';
            elseif nargin == 3
                model = 'other';
            end
            newData = TimeSeriesData_jy(filename, model);
            for i=1:length(tsd.wavelengthData)
                endTime = tsd.wavelengthData(i).times(end);
                % append the time array of part 2adding the final time of
                % part 1
                tsd.wavelengthData(i).times = [tsd.wavelengthData(i).times,...
                    (endTime + newData.wavelengthData(i).times)+gap/60];
                % append the data of part 2
                tsd.wavelengthData(i).data = [tsd.wavelengthData(i).data;...
                    newData.wavelengthData(i).data];
            end
            tsd.numberOfTimePoints = tsd.numberOfTimePoints + newData.numberOfTimePoints;
        end
        
        % set the name of samples in a column
        % 'lines' is an optional argument to set which lines to include
        % if 'lines' is not defined, then assume lines = 1:8
        function tsd = setSampleName(tsd, name, column, lines)
            if nargin == 3
                lines = 1:8;
            end
%             if ~isempty(find(tsd.name
            for i = lines,
                for j = 1:length(column)
                    well = [char(64 + i) num2str(column(j))];
                    tsd = setWell(tsd, name, well);
                end;
            end;
%             else
%               'Please use unique sample names.'
        end
        
        
        % assign a well to a given sample
        function tsd = setWell(tsd, sampleName, well)
            s = tsd.getSampleNumber(sampleName);
            % if sample doesn't exist yet, add it
            if (s == 0),
                tsd.samples(end+1).name  = sampleName;
                s = length(tsd.samples);
                tsd.samples(end).wells = [];
                tsd.samples(end).columns = [];
            end;
            tsd.samples(s).wells{end+1} = well;
            % convert wells to columns
            tsd.samples(s).columns(end+1) = tsd.findColumnNumber(well);
            % append data to the sample
            for w = 1:length(tsd.wavelengthData)
                tsd.samples(s).wavelength(w).data    =...
                    tsd.wavelengthData(w).data(:,...
                    tsd.samples(s).columns);
            end;
        end
        
        
        %%%%%%%%%% GET DATA from TimeSeriesData object
        
        % get a plate data for a given cycle number
        function data = getPlateMatrix(tsd, wavelength, cycle, iOrders)
            data = tsd.wavelengthData(wavelength).data(cycle, :);
            if nargin > 3
                data = data(iOrders);
            end
            data = reshape(data, 8, 12);
        end
        
        % get a plate data for a given cycle number
        % get only the center 6x10matrix. To be used with water around.
        function data = getPlateMatrixSmall(tsd, wavelength, cycle, iOrders)
            data = tsd.wavelengthData(wavelength).data(cycle, :);
            if nargin > 3
                data = data(iOrders);
            end
            data = reshape(data, 6, 10);
        end
        
        % get the matrix data for a specific wavelength and sample number
        function data = getData(tsd, wavelength, sampleNumber)
            data = tsd.samples(sampleNumber).wavelength(wavelength).data;
        end
        
        % get the array of times
        function times = getTimes(tsd, wavelength)
            times = tsd.wavelengthData(wavelength).times;
        end
        
        % get the data for a sample name
        function data = getDataFromName(tsd, wavelength, sampleName)
            % find the number of sample
            for s = 1:length(tsd.samples),
                if strcmp(tsd.samples(s).name, sampleName)
                    break
                end
            end
            % get the data
            data = tsd.getData(wavelength, s);
        end
        
        % search for the number of a sample with a given name
        function n = getSampleNumber(tsd, name)
            n = 0;
            for i = 1:length(tsd.samples)
                if strcmp(name, tsd.samples(i).name)
                    n = i;
                    return;
                end;
            end;
        end;
        
        % search for the column number that corresponds to a well
        function n = findColumnNumber(tsd, well)
            n = 0;
            for i = 1:length(tsd.map)
                if strcmp(well, tsd.map{i})
                    n = i;
                    return;
                end;
            end;
        end;
        
        
        %%%%%%%%%% Data correction functions
        
        % do blank correction using median over entire blak replicates
        function tsd = performBlankCorrection(tsd, blankName)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelengthData)
                % get the blankData
                blankData = tsd.getDataFromName(w, blankName);
                % calculate median along 2nd dimension (along replicates)
                blankData = median(blankData, 2);
                % subtract that from all data
                for s = 1:length(tsd.samples),
                    % get the data
                    data =  tsd.samples(s).wavelength(w).data;
                    % subtract blank to each column
                    for i = 1:size(data, 2)
                        data(:, i) = data(:, i) - blankData;
                    end
                    % rewrite the variable
                    tsd.samples(s).wavelength(w).data = data;
                end
            end
            tsd = tsd.updateSampleData;
        end
        
        % do blank correction but take as correction the same line in
        % the blank sample rather than median over entire blak replicates
        % This function has two ways to be used used.
        % If the argument 'samples' is not supplied, then the function
        % preforms blank correction for all the samples.
        % If 'samples' is supplied as a vector of sample numbers, then the
        % correction is carried out only for the samples in the array.
        function tsd = performBlankCorrectionPerLine(tsd, blankName,...
                samples)
            if (nargin == 2)
                samples = 1:length(tsd.samples);
            end;
            
            % determine if the blank is at the bottom or not
            if tsd.getSampleNumber(blankName) > 12 
                %Depends on the blank being declared LAST
                
            
                % If it is at the bottom, do as normal
                for w = 1:length(tsd.wavelengthData)
                    % get the blankData
                    blankData = tsd.getDataFromName(w, blankName);
                    % get median accross replicates
                    blankData = median(blankData, 2);
                    % subtract that from all data
                    for s = samples,
                        % get the data
                        data =  tsd.samples(s).wavelength(w).data;
                        for i = 1:size(data, 2)
                            data(:, i) =  data(:, i) - blankData;
                        end;
                        % rewrite the variable
                        tsd.samples(s).wavelength(w).data = data;
                    end
                end
            
            else %do the correction per row
                
                % do correction for all wavelengths
                for w = 1:length(tsd.wavelengthData)
                    % get the blankData
                    blankData = tsd.getDataFromName(w, blankName);
                    % get the number of blank sample
                    blankSample = tsd.getSampleNumber(blankName);
                    % subtract blanks from data, line by line
                    for s = samples,
                        % get the data
                        data =  tsd.samples(s).wavelength(w).data;
                        sampleWells = tsd.samples(s).wells;
                        for i = 1:size(data, 2)
                            sampleWell = sampleWells{i};
                            % find line in sample that corresponds to this
                            % line (first well name matches)
                            matching = strfind(tsd.samples(blankSample).wells,...
                                sampleWell(1));
                            match = not(cellfun(@isempty, matching));
                            % calculate median of blanks in the same row or column                        
                            data(:, i) =  data(:, i) -...
                                median(blankData(:, match), 2);
                        end;
                        % rewrite the variable
                        tsd.samples(s).wavelength(w).data = data;
                    end
                end
                
            end
            
            tsd = tsd.updateSampleData;
        end
        
        % correction using blank name
        function tsd=performDiffBlankCorrection(tsd, sampleindex, blankname)
            blankIdx = ismember({tsd.samples.name}, blankname);
            blankOD = tsd.samples(blankIdx).wavelength.data;
            blankODmedian = median(blankOD,2);
            for i=1:length(sampleindex)
                od = tsd.samples(sampleindex(i)).wavelength.data - ...
                                    repmat(blankODmedian, 1, size(tsd.samples(sampleindex(i)).wavelength.data, 2));
                tsd.samples(sampleindex(i)).wavelength.data = od;
            end
        end
        
        
        % correction is carried out only for the samples in the array.
        function tsd = subtractBackgroundValue(tsd, value,...
                samples, w1, w2)
            % subtract that from all data
            for s = samples,
                % get the data
                data =  tsd.samples(s).wavelength(w1).data;
                for i = 1:size(data, 2)
                    data(:, i) =  data(:, i) - value *...
                        tsd.samples(s).wavelength(w2).data(:, i);
                end;
                % rewrite the variable
                tsd.samples(s).wavelength(w1).data = data;
            end
            tsd = tsd.updateSampleData;
        end
        
        % do blank correction taking as background the first 10 point
        function tsd = performBlankCorrection10Points(tsd)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelengthData)
                % get the data
                data =  tsd.wavelengthData(w).data;
                % define the blankData
                blankData = median(data(1:10, :), 1);
                blankData = repmat(blankData, size(data, 1), 1);
                % subtract
                data =  data - blankData;
                % rewrite the variable
                tsd.wavelengthData(w).data = data;
            end
            tsd = tsd.updateSampleData;
        end
        
        % use after blanking to make sure sample data is up to date
        function tsd = updateSampleData(tsd)
            % subtract that from all data
            for s = 1:length(tsd.samples),
                for w = 1:length(tsd.wavelengthData)
                    tsd.wavelengthData(w).data(:,...
                        tsd.samples(s).columns)=...
                        tsd.samples(s).wavelength(w).data;
                end
            end
        end
        
        % move the curves in the x-axis by a time step tau
        function tsd = moveInTime(tsd, samples, tau)
            % do correction for all wavelengths
            for w = 1:length(tsd.wavelengthData)
                time = tsd.getTimes(w);
%                 indexValid = find( and(time>tau, time<(time(end)-tau)) );
                indexValid = find(time>tau);
                % correct all samples in array
                for s = samples,
                    % get the data
                    data =  tsd.samples(s).wavelength(w).data;
                    for i = 1:size(data, 2)
                        newTime  = time(indexValid)-tau;
                        dataMoved = interp1(newTime, data(indexValid, i), time);
                        data(:, i) =  dataMoved;
                    end;
                    % rewrite the variable
                    tsd.samples(s).wavelength(w).data = data;
                end
            end
        end;
        
        
        %%%%%%%%%% Plotting functions
        
        
        % plot the time series of the median for a given sample
        function h = plotMedian(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            h(1) = plot(times, median(data, 2), 'b-');
            set(h(1), 'Linewidth', 3);
        end
        
        % plot the time series of ranges for a given sample
        function h = plotRangesAsLines(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            
%             if wavelength == length(tsd.wavelengthData);%Assumes the last 
                %wavelength is the one needing autocorrection
            if findstr(lower(tsd.wavelengthData(wavelength).name), 'gfp')
                [rangeU, rangeL] = tsd.calculateGFPRanges(sampleNumber);
                
                h(1) = plot(times, rangeU, 'k-');
                hold on;
                h(2) = plot(times, rangeL, 'k-');
                hold off;
               
            else      
                h(1) = plot(times, max(data, [], 2), 'k-');
                hold on;
                h(2) = plot(times, min(data, [], 2), 'k-');
                hold off;
            end
            
        end
        
        function [rangeU, rangeL] = calculateGFPRanges(tsd,sampleNumber)
        
            if size(tsd.samples) > 4 % If the sample may have auto correction present
                
                if ~isempty(tsd.samples.af)
                    gfpData = tsd.getData(2, sampleNumber);
                    % subtract that from all data
                    for i = 1:size(gfpData, 2)
                        gfpDataU(:, i) =  gfpData(:, i) - myAutoData.samples(sampleNumber).af(:,1);
                        gfpDataL(:,i) = gfpData(:, i) - myAutoData.samples(sampleNumber).af(:,2);
                    end;
                    
                    rangeU = max(gfpDataU,[],2);
                    rangeL = min(gfpDataL,[],2);
                    
                else
                    gfpData = tsd.getData(2, sampleNumber);
                    rangeU = max(gfpData, [], 2);
                    rangeL = min(gfpData, [], 2);
                    
                end
                
            else % No possible autocorrection data present
                gfpData = tsd.getData(2, sampleNumber);
                rangeU = max(gfpData, [], 2);
                rangeL = min(gfpData, [], 2);
                
            end
                    
        end
        
        % plot the time series of ranges for a given sample
        function h = plotRangesAsArea(tsd, wavelength, sampleNumber)
            times = tsd.getTimes(wavelength);
            data  = tsd.getData(wavelength, sampleNumber);
            
            if wavelength ~= length(tsd.wavelengthData);%Assumes the last 
                %wavelength is the one needing autocorrection
                maxData = max(data, [], 2);
                minData = min(data, [], 2);
            else
                [maxData, minData] = tsd.calculateGFPRanges(sampleNumber);
            end
             
            range   = maxData - minData;
            % take out NaN's, otherwise area won't work properly
            nans = isnan(minData) | isnan(range) | (minData <= 0);
            minData = minData(~nans);
            range = range(~nans);
            times = times(~nans);
            % plot the area
            h = area(times, [minData, range]);
            set(h, 'EdgeColor', 'none');
            % erase the first area
            set(h(1), 'FaceColor', 'none');
            % return only the handle to the second area
            h = h(2);
        end
        
        
        % plot the time series of the median for a given sample
        % with ranges overlapped
        function h = plotMedianWithRanges(tsd, wavelength, sampleNumber)
%             times = tsd.getTimes(wavelength);
%             data  = tsd.getData(wavelength, sampleNumber);
            h(1)= tsd.plotMedian(wavelength, sampleNumber);
            hold on;
            h2 = tsd.plotRangesAsLines(wavelength, sampleNumber);
            hold off;
            h = [h h2(1) h2(2)];
            set(h(1), 'Linewidth', 3);
            set(gca, 'YScale', 'log');
            title(tsd.samples(sampleNumber).name);
            set(gca,  'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelengthData(wavelength).name);
        end
        
        
        % plots OD and GFP in the same plot, normalized by the entire
        % sample (assumes OD is wavelength 1 and GFP is wavelength 2)
        function plotNormalizedODAndGFP(tsd, sampleNumber)
            %
            odData = (tsd.getData(1, sampleNumber));
            medianOd = median(odData, 2);
            maxVal = max(medianOd(~isnan(medianOd)));
            minVal = 0.01;
            odData = (odData - minVal) ./ (maxVal - minVal);
            
            % gfp (second wavelength)
            gfpData = tsd.getData(2, sampleNumber);
            medianGfp = median(gfpData, 2);
            maxVal = max(medianGfp(~isnan(medianGfp)));
            minVal = 5;
            gfpData = (gfpData - minVal) ./ (maxVal - minVal);
            
            % time
            times = tsd.getTimes(1);
            
            % Hilary
            startVal = length(times)-length(odData)+1; %The +1 is to compensate for indexing beginning at 1
            timesPlot = times(startVal:end);
            
            maxTime = max(timesPlot(~isnan(timesPlot)));
            minTime = min(timesPlot(~isnan(timesPlot)));
            plot(timesPlot, median(odData, 2), 'k', 'LineWidth', 2);
            hold on;
            plot(timesPlot, median(gfpData, 2),...
                'Color', [0 0.8 0], 'LineWidth', 2);
            plot(timesPlot, max(odData, [], 2), 'k-');
            plot(timesPlot, min(odData, [], 2), 'k-');
            plot(timesPlot, max(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(timesPlot, min(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            hold off;
            set(gca, 'YLim', [0 1],...
                'XLim', [minTime maxTime], 'YTickLabel', []);
        end
        
        % plots GFP data from multiple samples in the same plot, 
        % normalized by the entire sample (GFP is wavelength 2)
        function plotNormalizedGFP(tsd, sampleNumbers)
            % time
            times = tsd.getTimes(1);
            numWells = min(size(tsd.getData(2,sampleNumbers(1))));
            
            %pre-allocate array for normalized gfp data
            allGFP = zeros(length(times), length(sampleNumbers), numWells);
            
            % get all gfp data
            for i = 1:length(sampleNumbers)
                gfpData = tsd.getData(2, sampleNumbers(i));
                medianGfp = median(gfpData, 2);
                maxVal = max(medianGfp(~isnan(medianGfp)));
                minVal = 5;
                gfpData = (gfpData - minVal) ./ (maxVal - minVal);
%                 max(gfpData)
%                 min(gfpData)
                allGFP(:,i,:) = gfpData;
                clear gfpData
            end
            
            startVal = length(times)-length(allGFP(:,1))+1; %The +1 is to compensate for indexing beginning at 1
            timesPlot = times(startVal:end);
            
            maxTime = max(timesPlot(~isnan(timesPlot)));
            minTime = min(timesPlot(~isnan(timesPlot)));
            
            cmap = jet(length(sampleNumbers));
            
            for j = 1:length(sampleNumbers)
                plot(timesPlot, median(allGFP(:,j,:), 3), 'Color', cmap(j,:), 'LineWidth', 2);
                hold on;
%                 plot(timesPlot, max(max(allGFP(:,j), [], 2)), '-', 'Color', cmap(j,:));
%                 plot(timesPlot, min(allGFP(:,j), [], 2), '-', 'Color', cmap(j,:));
                plot(timesPlot, max(allGFP(:,j,:), [], 3), '-', 'Color', cmap(j,:));
                plot(timesPlot, min(allGFP(:,j,:), [], 3), '-', 'Color', cmap(j,:));
            end
            hold off;
            set(gca, 'YLim', [0 1],...
                'XLim', [minTime maxTime]);%, 'YTickLabel', []);
        end
        
        
        % plots three wavelengths in the same plot, normalized by the entire
        % sample. the example here is for od, gfp and pyoverdine
        % (assumes OD is wavelength 1, GFP is wavelength 2 and pyoverdine is 3)
        function plotThreeWavelengths(tsd, sampleNumber)
            % od (first wavelength)
            odData = (tsd.getData(1, sampleNumber));
            maxVal = max(odData(~isnan(odData)));
            minVal = min(odData(~isnan(odData)));
            odData = (odData - minVal) ./ (maxVal - minVal);
            % gfp (second wavelength)
            gfpData = tsd.getData(2, sampleNumber);
            maxVal = max(gfpData(~isnan(gfpData)));
            minVal = min(gfpData(~isnan(gfpData)));
            gfpData = (gfpData - minVal) ./ (maxVal - minVal);
            % pyo (third wavelength)
            pyoData = tsd.getData(3, sampleNumber);
            maxVal = max(pyoData(~isnan(pyoData)));
            minVal = min(pyoData(~isnan(pyoData)));
            pyoData = (pyoData - minVal) ./ (maxVal - minVal);
            % time
            times = tsd.getTimes(1);
            maxTime = max(times(~isnan(times)));
            minTime = min(times(~isnan(times)));
            plot(times, median(odData, 2), 'k', 'LineWidth', 2);
            hold on;
            plot(tsd.getTimes(2), median(gfpData, 2),...
                'Color', [0 0.8 0], 'LineWidth', 2);
            plot(tsd.getTimes(3), median(pyoData, 2),...
                'Color', [0 0 1], 'LineWidth', 2);
            plot(times, max(odData, [], 2), 'k-');
            plot(times, min(odData, [], 2), 'k-');
            plot(times, max(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(times, min(gfpData, [], 2), '-', 'Color', [0 0.8 0]);
            plot(times, max(pyoData, [], 2), '-', 'Color', [0 0 1]);
            plot(times, min(pyoData, [], 2), '-', 'Color', [0 0 1]);
            hold off;
            set(gca, 'YLim', [0 1],...
                'XLim', [minTime maxTime], 'YTickLabel', []);
        end
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        % last argument (rangeFlag) is optional:
        %  rangeFlag = true (default) plots ranges as thin lines
        %  rangeFlag = false plots only the median
        function plotManySamples(tsd, wavelength, sampleNumbers, rangeFlag)
            % construct a colormap
            labels = [];
            cmap = jet(length(sampleNumbers));
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotMedian(wavelength, sampleNumbers(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(sampleNumbers(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
            if ((nargin == 4) && (rangeFlag ~= false)) || (nargin == 3)
                for i = 1:length(sampleNumbers)
                    hold on;
                    %h = tsd.plotRangesAsArea(wavelength, sampleNumbers(i));
                    %set(h, 'FaceColor', cmap(i, :));
                    h = tsd.plotRangesAsLines(wavelength, sampleNumbers(i));
                    set(h, 'Color', cmap(i, :));
                end;
            end
            set(gca, 'YScale', 'log', 'XLim', [0 24], 'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelengthData(wavelength).name);
        end;
        
        % plot wavelength w2 over wavelength w2 Medians
        % for sample s Kerry's code
        function h = plotW1OverW2Median(tsd, w1, w2, s)
            
            times = tsd.getTimes(w1);
            d1 = tsd.getData(w1, s);
            d2 = tsd.getData(w2, s);
            Ov = d1./d2;
            MedOv = median(Ov, 2);
            wavelength = Ov;
            h = plot(times(1,1:length(d1)), MedOv, '-', 'LineWidth', 2);
            hold on
           
            xlabel('Time[h]');
            ylabel([tsd.wavelengthData(w1).name '/'...
                tsd.wavelengthData(w2).name]);
         
           
        end;
        
        % same as plotW1OverW2Median but works with array
        % of samples Kerry's code
        function plotW1OverW2MedianMany(tsd, w1, w2, samples)
            
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            
            for i = 1:length(samples)
                hold on;
                h = tsd.plotW1OverW2Median(w1, w2, samples(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            
            legend(labels, 'Location', 'Best');
         
        end
        
        
        %Plot rhlAB promoter activity
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Kerry's Code 09/02/2014%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = plotPromoter(tsd, w1, w2, s)
            d1 = tsd.getData(w1, s);
            d2 = tsd.getData(w2, s);
            times = tsd.getTimes(w1);
            % smooth the data line
             odfilter = TimeSeriesData_jy.filterData2(median(d1,2),3);
             gfpfilter = TimeSeriesData_jy.filterData2(median(d2,2),3);

             dgfp = tsd.calculateRate(times, gfpfilter);
             dgfp = dgfp./odfilter;
                     
          v1 = find(odfilter > 0.01);                   
          h = plot(times(v1), dgfp(v1), 'k-', 'Linewidth', 3);
         
         % h = plot(times, dgfp .* v1, '-', 'Linewidth', 2);

            xlabel('Time [hours]');
            ylabel('activity');    
    
        end
        
        
        function plotPromoterMany(tsd, w1, w2, samples)
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            
            for i = 1:length(samples)
                hold on;
                h = tsd.plotPromoter(w1, w2, samples(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            legend(labels, 'Location', 'Best');
        end
        
        
        
        % plot wavelength w2 as function of wavelength w2
        % for sample s
        function h = plotWavelength1VsWavelength2(tsd, w1, w2, s)
            d1 = tsd.getData(w1, s);
            d2 = tsd.getData(w2, s);
            h = plot(d1(:), d2(:), '+');
            xlabel(tsd.wavelengthData(w1).name);
            ylabel(tsd.wavelengthData(w2).name);
        end;
        
        % same as plotWavelength1VsWavelength2 but works with array
        % of samples
        function plotW1VsW2ManySamples(tsd, w1, w2, samples)
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            for i = 1:length(samples)
                hold on;
                h = tsd.plotWavelength1VsWavelength2(w1, w2, samples(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            legend(labels, 'Location', 'Best');
        end
        
        % plts the specific growth rate time series
        function [d1, h] = plotSpecificGrowthW1(tsd, w1, s)
            d1 = tsd.getData(w1, s);
            time = tsd.getTimes(w1); %***
            d1 = TimeSeriesData.calculateSpecificRate(time(1:length(d1)), d1);
            h = plot(time, median(d1, 2), '-');
            ylabel(['growth' tsd.wavelengthData(w1).name]);
            xlabel('time [h]');
        end;
        
        function h = plotSpecificGrowthW1VsW2(tsd, w1, w2, s, thresholdW1, useMedian)
            % default useMedian = 0
            if nargin == 5
                useMedian = 0;
            end
            d1 = tsd.getData(w1, s);
            % take only the median of w1
            if useMedian == 1
                d1 = median(d1, 2);
            end;
            %
            time = tsd.getTimes(w1);
            d1 = TimeSeriesData.calculateSpecificRate(time, d1);
            d2 = TimeSeriesData.filterData(tsd.getData(w2, s));
            % take only the median of w2
            if useMedian == 1
                d2 = median(d2, 2);
            end
            %
            % threshold for wvelength w1
            d1Filterd = TimeSeriesData.filterData(tsd.getData(w1, s));
            % take only the median of w1
            if useMedian == 1
                d1Filterd = median(d1Filterd, 2);
            end
            %
            d1(d1Filterd < thresholdW1) = [];
            d2(d1Filterd < thresholdW1) = [];
            d1Filterd(d1Filterd < thresholdW1) = [];
            h = plot(d1(:), d2(:)./d1Filterd(:), '.');
            ylabel(tsd.wavelengthData(w2).name);
            % use second line to normalize gfp/od600
            %h = plot(d1(:), d2(:)./d1Filterd(:), '+');
            %ylabel([tsd.wavelengthData(w2).name '/'...
            %    tsd.wavelengthData(w1).name]);
            xlabel(['specific growth ' tsd.wavelengthData(w1).name]);
        end;
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        function plotSpecificGrowthW1VsW2Many(tsd, w1, w2, samples, thresholdW1, useMedian)
            % default useMedian = 0
            if nargin == 5
                useMedian = 0;
            end
            
            % construct a colormap
            labels = [];
            cmap = jet(length(samples));
            for i = 1:length(samples)
                hold on;
                h = tsd.plotSpecificGrowthW1VsW2(w1, w2, samples(i), thresholdW1, useMedian);
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(samples(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
        end
        
        % plot sample numbers provided in a vector with ranges
        % and add legend
        function plotManySamplesWithArea(tsd, wavelength, sampleNumbers)
            % construct a colormap
            labels = [];
            cmap = jet(length(sampleNumbers));
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotMedian(wavelength, sampleNumbers(i));
                set(h, 'Color', cmap(i, :));
                labels{end+1} = tsd.samples(sampleNumbers(i)).name;
            end;
            legend(labels, 'Location', 'SouthEast');
            for i = 1:length(sampleNumbers)
                hold on;
                h = tsd.plotRangesAsArea(wavelength, sampleNumbers(i));
                set(h, 'FaceColor', cmap(i, :));
            end;
            set(gca, 'YScale', 'log', 'XLim', [0 24], 'YLim', [0.01 1]);
            xlabel('time [h]');
            ylabel(tsd.wavelengthData(wavelength).name);
        end;
        
        % plot a matrix representing data from a plate
        % use flagForSmall if plate has water all around
        function plotMatrix(tsd, wavelength, cycleNumber, column, line,...
                flagForSurface, iOrders, flagForSmall)
            if nargin > 6
                if flagForSmall
                    matrix2Plot =...
                        tsd.getPlateMatrixSmall(wavelength, cycleNumber, iOrders);
                else
                    matrix2Plot =...
                        tsd.getPlateMatrix(wavelength, cycleNumber, iOrders);
                end
            else
                matrix2Plot = tsd.getPlateMatrix(wavelength, cycleNumber);
            end
            if flagForSurface
                surf(matrix2Plot);
            else
                imagesc(matrix2Plot);
            end
            times = tsd.getTimes(wavelength);
            title(sprintf('%s at %0.1f h',...
                tsd.wavelengthData(wavelength).name,...
                times(cycleNumber)));
            colorbar;
            set(gca, 'XTick', 1:size(matrix2Plot, 2),...
                'XTickLabel', line);
            set(gca, 'YTick', 1:size(matrix2Plot, 1),...
                'YTickLabel', column);
        end;
        
        %%%%%%%%%% Growth curve synchronization functions
        
        % find the array of time delays between sample 'referenceSample'
        % and the samples in array sampleNumbers
        function tauArray = computeTimeDelays(tsd,...
                wavelength, referenceSample, sampleNumbers)
            time = tsd.getTimes(wavelength);
            tauArray = zeros(1, length(sampleNumbers));
            curve1 = tsd.getData(wavelength, referenceSample);
            for i = 1:length(sampleNumbers)
                s = sampleNumbers(i);
                curve2 = tsd.getData(wavelength, s);
                fToMin =...
                    @(tau)...
                    (TimeSeriesData.alignmentError(time,...
                    curve1, curve2, tau));
                tauArray(i) = fminbnd(fToMin, 0, 24,...
                    optimset('TolX', 1e-2));
            end
        end
        
        % find the array of time delays between sample 'referenceSample'
        % and the samples in array sampleNumbers
        function tauMatrix = computeTimeDelayMatrix(tsd,...
                wavelength, sampleNumbers)
            n = length(sampleNumbers);
            tauMatrix = zeros(n, n);
            for i = 1:n
                for j = (i+1):n
                    tauMatrix(i, j) = tsd.computeTimeDelays(...
                        wavelength, sampleNumbers(i), sampleNumbers(j));
                end
            end
        end
        
        % calulate a total error of overlap for a given array of delays
        function errorVal = computeErrorMatrix(tsd,...
                wavelength, sampleNumbers, tauArray)
            tauArray = [0 tauArray];
            n = length(sampleNumbers);
            errorMatrix = zeros(n, n);
            time = tsd.getTimes(wavelength);
            for i = 1:n
                for j = (i+1):n
                    curve1 = tsd.getData(wavelength, sampleNumbers(i));
                    curve2 = tsd.getData(wavelength, sampleNumbers(j));
                    tau = tauArray(j) - tauArray(i);
                    errorMatrix(i, j) =...
                        TimeSeriesData.alignmentError(time, curve1, curve2, tau);
                end
            end
            errorVal = sum(errorMatrix(:));
        end
        
        
        % calculate an array of delays by minimizing the
        % total error of overlap between each pair of curves
        function tauArray = optimizeTauArray(tsd,...
                wavelength, sampleNumbers)
            tauArray = tsd.computeTimeDelays(wavelength,...
                sampleNumbers(1), sampleNumbers(2:end));
            fToMin =...
                @(x)...
                (computeErrorMatrix(tsd,...
                wavelength, sampleNumbers, x));
            tauArray = fminsearch(fToMin, tauArray);
            
        end
        
        
        %%%%%%%%%% Growth phase analysis
        function tsd = setPhaseTimes(tsd, samples, phaseTimes)
            if (length(samples) ~= size(phaseTimes, 2))
                error(...
                    'length of "samples" must match number of columns in "phaseTimes"');
            end
            % initialize phaseTimes
            if isempty(tsd.phaseTimes)
                tsd.phaseTimes = zeros(3, length(tsd.samples));
            end
            % set the phase time
            for i = 1:length(samples)
                tsd.phaseTimes(:, samples(i)) = phaseTimes(:, i);
            end;
        end
        
        
        
        function plotGrowthPhase (tsd, wavelength, samples, timeRangeMax, logFlag)
            % Construct a colormap for the 3 phases
            cmap = jet(length(samples)); %3 phases * 3 strains - DONT COMPARE WITHIN STRAINS 8/1/13 TO DO - FIX THIS
            % strainNum = ceil(sampleNumber/4);
            
            % Setup figure to plot d(wavelength data)/dt, where dt = 10 min
            % = 1/6 hrs
            % Use phaseTimes (a 1x3 matrix) to determine where to split the
            % wavelength data into 4 segments (3 phases + lag phase). Plot each of the 3 segments
            
            % Hardcoding phaseTimes if needed
            % First element is where phase 1 starts, second element is where phase 2 starts
            % phaseTimes = [14.5 21.75 26.25]; %8/1/13 - based on Hilary's growth curve with glucose and cbra mutant 7/29/13
            % phaseTimes = [12, 18];
            xshift = linspace(-0.3, 0.3, length(samples));
            if length(samples) > 1
                deltax = xshift(2) - xshift(1);
            else
                deltax = 0.1;
            end
            for i = 1:length(samples)
                phaseTimes = tsd.phaseTimes(:, samples(i));
                
                timeData = tsd.getTimes(wavelength);
                sampleData = tsd.getData(wavelength, samples(i));
                timeData = timeData(1:length(sampleData));
                deltaT = (timeData(2)-timeData(1)); %This is a slight approximation
                medianData = median(sampleData,2);
                filteredMedianData = tsd.filterData(medianData);
                logData = log(filteredMedianData);
                                
                %preallocation of phaseMatrix
                phase = zeros(1,length(timeData));
                
                phase(timeData < timeRangeMax) = 3;
                phase(timeData < phaseTimes(3)) = 2;
                phase(timeData < phaseTimes(2)) = 1;
                phase(timeData < phaseTimes(1)) = 0;
                
                
                %calculate the change in gfp               
                if logFlag == true
                    filteredMedianDData = [0; diff(logData)/deltaT];
                else
                    filteredMedianDData =...
                        [0; diff(filteredMedianData)/deltaT];
                end
                
                phaseR = phase + xshift(i) + rand(size(phase))*deltax/2;
                
                scatter(phaseR(phase>0)', filteredMedianDData(phase>0)', 40,...
                    cmap(i, :), 'LineWidth', 1.5)% 'MarkerEdgeColor', c)
                hold on;
            end
            hold off;
            xlabel('Phase')
            ylabel(strcat('d(',char(tsd.wavelengthData(wavelength).name),')/dt'));
        end
        
        %%%%%%%%%% Autofluorescence correction
        
        
        
        % do blank correction but take as correction the same line in
        % the blank sample rather than median over entire blak replicates
        % This function has two ways to be used used.
        % If the argument 'samples' is not supplied, then the function
        % preforms blank correction for all the samples.
        % If 'samples' is supplied as a vector of sample numbers, then the
        % correction is carried out only for the samples in the array.
        function tsd = performAutoCorrection(tsd, autoSample, sample, w)
            % get the blankData
            autoData = tsd.getData(w, autoSample);
            % get median accross replicates
            autoData = median(autoData, 2);
            % subtract that from all data
            % get the data
            data =  tsd.samples(sample).wavelength(w).data;
            for i = 1:size(data, 2)
                data(:, i) =  data(:, i) - autoData;
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(w).data = data;
            tsd = tsd.updateSampleData;
        end
        
        function tsd = performAuto380(tsd, sample, od380w, gfpw)
            % get the od380 data
            od380Rep = median(tsd.getData(od380w, sample),2);
            auto380 = 10^3*[1.0480 0.0491]; %robust for the 12092013 dataset
            % create the auto correction data (median only)
            myAutoData380 = od380Rep*auto380(1)+auto380(2);
            % get the gfp data
            gfpData = tsd.getData(gfpw, sample);
            % subtract that from all data
            for i = 1:size(gfpData, 2)
                gfpData(:, i) =  gfpData(:, i) - myAutoData380;
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(gfpw).data = gfpData;
            tsd = tsd.updateSampleData;
        end
        
        function tsd = performAuto380Log(tsd, sample, od380w, gfpw)
            % get the od380 data
            od380Rep = median(tsd.getData(od380w, sample),2);
            auto380 = [6.9060, 0.7173];%[7.0914 0.7864]; % CHECK THESE
            % create the auto correction data (median only)
            myAutoData380 = exp(auto380(1) + log(od380Rep) * auto380(2));
            % get the gfp data
            gfpData = tsd.getData(gfpw, sample);
            % subtract that from all data
            for i = 1:size(gfpData, 2)
                gfpData(:, i) =  gfpData(:, i) - myAutoData380;
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(gfpw).data = gfpData;
            tsd = tsd.updateSampleData;
        end
        
        function tsd = performAuto600Log(tsd, sample, od600w, gfpw)
            % get the od600 data
%             tsd.sample(sample).addError(:,1) = [];
%             tsd.sample(sample).addError(:,2) = [];
            od600Rep = median(tsd.getData(od600w, sample),2);
            auto600 = [8.2034, 0.9734];
            auto600L = [8.1872, 0.9672];
            auto600U = [8.2197, 0.9796];
            
            tsd.samples(sample).af = [];
            
            % create the auto correction data (median only)
            myAutoData600 = exp(auto600(1) + log(od600Rep) * auto600(2));
            myAutoData600U = exp(auto600U(1) + log(od600Rep) * auto600U(2));
            myAutoData600L = exp(auto600L(1) + log(od600Rep) * auto600L(2));
            
            tsd.samples(sample).af = [myAutoData600U, myAutoData600L];

            % get the gfp data
            gfpData = tsd.getData(gfpw, sample);
            % subtract that from all data
            for i = 1:size(gfpData, 2)
                gfpData(:, i) =  gfpData(:, i) - myAutoData600;
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(gfpw).data = gfpData;
            tsd = tsd.updateSampleData;
        end
        
        function tsd = performAuto600LogW(tsd, sample, od600w, gfpw)

            
            auto600 = [8.2034, 0.9734];
%             auto600L = [8.1872, 0.9672];
%             auto600U = [8.2197, 0.9796];
            
%             tsd.samples(sample).af = [];
            
            % get the od600 data
            od600Rep = tsd.getData(od600w, sample);
            logoddata = log(od600Rep);
            
            % get the sample wells
            wells = tsd.samples(sample).wells;

            % get the gfp data
            gfpData = tsd.getData(gfpw, sample);
            
            %preallocation
            myAuto600 = zeros(length(tsd.getTimes(1)));
            % subtract AF from each well in the sample
            for i = 1:length(wells)
                myAuto600(:, i) = exp(auto600(1) + logoddata(:, i) * auto600(2));
                gfpData(:, i) =  gfpData(:, i) - myAuto600(:,i);
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(gfpw).data = gfpData;
            tsd = tsd.updateSampleData;
        end
        
        function tsd = performAuto600LogGlu(tsd, sample, od600w, gfpw)
            % get the od600 data
%             tsd.sample(sample).addError(:,1) = [];
%             tsd.sample(sample).addError(:,2) = [];
            od600Rep = median(tsd.getData(od600w, sample),2);
            auto600 = [7.4255, 0.9602];
            auto600L = [7.3391, 0.9302];
            auto600U = [7.5120, 0.9901];
            
            tsd.samples(sample).af = [];
            
            % create the auto correction data (median only)
            myAutoData600 = exp(auto600(1) + log(od600Rep) * auto600(2));
            myAutoData600U = exp(auto600U(1) + log(od600Rep) * auto600U(2));
            myAutoData600L = exp(auto600L(1) + log(od600Rep) * auto600L(2));
            
            tsd.samples(sample).af = [myAutoData600U, myAutoData600L];

            % get the gfp data
            gfpData = tsd.getData(gfpw, sample);
            % subtract that from all data
            for i = 1:size(gfpData, 2)
                gfpData(:, i) =  gfpData(:, i) - myAutoData600;
            end;
            % rewrite the variable
            tsd.samples(sample).wavelength(gfpw).data = gfpData;
            tsd = tsd.updateSampleData;
        end
        
        % runs performAutoCorrection for many samples
        % length of "samples"must be multiple of length of "autoSamples"
        function tsd = performAutoCorrectionManySamples(tsd, autoSamples, samples, w)
            if mod(length(samples), length(autoSamples)) ~= 0
                error('length of "samples" must be multiple of length of "autoSamples"');
            end
            for i = 1:length(samples)
                tsd = tsd.performAutoCorrection(autoSamples(mod(i-1, length(autoSamples))+1),...
                    samples(i), w);
            end
        end
        
        function tsd = performAuto380ManySamples(tsd, samples, od380w, gfpw)
            if nargin < 5
                gfpw = 3;
                if nargin < 4
                    od380w = 2;
                end
            end
            
            for i = 1:length(samples)
                tsd = tsd.performAuto380(samples(i), od380w, gfpw);
            end
        end
        
        function tsd = performAuto380ManySamplesLog(tsd, samples, od380w, gfpw)
            if nargin < 5
                gfpw = 3;
                if nargin < 4
                    od380w = 2;
                end
            end
            
            for i = 1:length(samples)
                tsd = tsd.performAuto380Log(samples(i), od380w, gfpw);
            end
        end
        
        function tsd = performAuto600ManySamplesLogW(tsd, samples, od600w, gfpw)
            if nargin < 4
                gfpw = length(tsd.samples(samples(1)).wavelength);
                if nargin < 3
                    od600w = 1;
                end
            end
           
            for i = 1:length(samples)
                tsd = tsd.performAuto600LogW(samples(i), od600w, gfpw);
            end
            
        end
        
        function tsd = performAuto600ManySamplesLog(tsd, samples, od600w, gfpw)
            if nargin < 4
                gfpw = length(tsd.samples(samples(1)).wavelength);
                if nargin < 3
                    od600w = 1;
                end
            end
           
            for i = 1:length(samples)
                tsd = tsd.performAuto600Log(samples(i), od600w, gfpw);
            end
            
        end
        function tsd = performAuto600ManySamplesLogGlu(tsd, samples, od600w, gfpw)
            if nargin < 4
                gfpw = length(tsd.samples(samples(1)).wavelength);
                if nargin < 3
                    od600w = 1;
                end
            end
           
            for i = 1:length(samples)
                tsd = tsd.performAuto600LogGlu(samples(i), od600w, gfpw);
            end
            
        end

        %%%%%%%%%% Lag correction
        
        
        function setLagTimes(tsd, lagTimes)
            tsd.lagTimes = lagTimes;
        end
        
        function lagTimes = getLagTimes(tsd)
            lagTimes = tsd.lagTimes;
        end
        
        
        % correct lag phase
        %Step 1 - find the first time each of the median samples comes
        %above a threshold .01
        %Step 2 - move all curves in time to the earliest time the median
        %samples go above their respective threshold, i.e. compute the
        %lag between each of the samples and do a lag shift.
        function tsd = autoLag(tsd, sampleNumbers, threshold)
            if nargin < 3
                threshold = 0.01; %this works for both od and gfp since gfp is usually very high
            end
            % calcute lag
            lagArray = tsd.computeLag(1, sampleNumbers, threshold);
            % shift data
            for i = 1:length(sampleNumbers)
                tsd = tsd.moveInTime(sampleNumbers(i), lagArray(2, i));
            end
        end

        
        function lag = computeLag(tsd,wavelength,sampleNumbers, threshold)
            % structure of output matrix lag =
            % [sample_1 od, sample2 od, ..., sample_n od]
            % [sample_1 t , sample2 t , ..., sample_n t ]
            % [sample_1 t-index, ..., sample_n t-index]
            
            %Preallocation and getting the time data - bases preallocation
            %size on the size of the data of the first dataset.
            times = tsd.getTimes(wavelength);
            lag = zeros(3,length(sampleNumbers));
            
            for i = 1:length(sampleNumbers)
                medianData = median(tsd.getData(1,sampleNumbers(i)),2);
                medianData(1) = 0; 
                validPoints = find(medianData > threshold);   
                
                lag(1,i) = medianData(validPoints(1));
                lag(2,i) = times(validPoints(1));
                lag(3,i) = validPoints(1); %track the index

            end
        end
        
        function tsd = forceGrowthCurveToPassInPoint(tsd, wavelength,...
                sampleNumbers, timeToCross, valueToCross)
            t  = tsd.getTimes(wavelength)';            
            for i = 1:length(sampleNumbers)
                s = sampleNumbers(i);
                medianData = median(tsd.getData(1,s),2);
                medianData(1) = 0; % take out first point
            
                firstPoint = find(medianData<(valueToCross/2), 1, 'last' ) +1 ;
                lastPoint = find(medianData>(valueToCross*2), 1 ) -1 ;
                v = firstPoint:lastPoint;
                try
                    mdl = LinearModel.fit(t(v),log(medianData(v)),'linear',...
                        'RobustOpts','on');
                    b = double(mdl.Coefficients(:,1));
                    deltaT = (log(valueToCross)-b(1))/b(2) - timeToCross;
                    tsd = tsd.moveInTime(s, deltaT);
                catch
                    warning(['regression did not work for sample '...
                        tsd.samples(s).name]);
                end
            end
        end
        
        
        % Designate the lagged Data for a sample number
        function tsd = setData(tsd, wavelength, sampleNumber, data)
            tsd.samples(sampleNumber).wavelength(wavelength).data = [];
            tsd.samples(sampleNumber).wavelength(wavelength).data = data;
        end
        
        
        %%% Jinyuan's code
        % get the matrix of GFP to od, return a 3D matrix
        function GfpOd = getGfpToOdmany(tsd, s) 
            
            for i=s
                Gfp = tsd.getData(2, i);
                Gfp = TimeSeriesData_jy.filterData(Gfp);
                Gr = tsd.getData(1, i);
                Gr = TimeSeriesData_jy.filterData(Gr);
                GfpOd(:,:,i) = Gfp ./ Gr;
                    % check to see if the gfp to od ratio is noisy
            %     time = tsd3.wavelengthData(1).times;
            %     figure
            %     plot(time, GfpOd / 1e4);
            %     title([tsd.samples(i).name ' GFP/OD (10^4)'], 'fontsize', 16)
            %     set(gca, 'ylim', [-1 5])
            end
        end
        
        % plot gfpod ratio for many samples
        function plotManyGfpOd(tsd, GfpOd, s)
            time = tsd.wavelengthData(1).times;
            mycolormap = jet(length(s));
            le = {};
            for i = 1:length(s)
                plot(time, median(GfpOd(:,:,s(i)),2) / 1e3, 'color', mycolormap(i,:), 'linewidth', 3);
                hold on
                le{end+1} = tsd.samples(s(i)).name;
            end
            legend(le)
            for i = 1:length(s)
                plot(time, max(GfpOd(:,:,s(i)),[],2) / 1e3, 'color', mycolormap(i,:), 'linewidth', 1);
                hold on
                plot(time, min(GfpOd(:,:,s(i)),[],2) / 1e3, 'color', mycolormap(i,:), 'linewidth', 1);
                hold on
            end
            hold off
            set(gca, 'ylim', [0 6], 'xlim', [0 24], 'fontsize', 14);
            title('Gfp/OD', 'fontSize', 16);
        end
    end % methods
    
    
    methods(Static)
        % extract the text in cell 'i' in a line of the spreadsheet
        function t = textInCell(line, i)
            commas = strfind(line, ','); % find indexes of commas
            if i == 1
                if isempty(commas)
                    t = line;
                else
                    t = line(1:commas(i)-1);
                end
            else
                t = line(commas(i-1)+1:commas(i)-1);
            end
        end % end textInCell
        
        
        % extract the number in cell 'i' in a line of the spreadsheet
        function a = checkIfLineHasData(line)
            commas = strfind(line, ','); % find indexes of commas
            a = ~isempty(commas);
            % check if there is time data at entry 2
            try
                if isnan(TimeSeriesData.numberInCell(line, 2))
                    a = false;
                end
            catch
                a = false;
            end
        end % end textInCell
        
        
        % extract the number in cell 'i' in a line of the spreadsheet
        function n = numberInCell(line, i)
            n = str2double(TimeSeriesData.textInCell(line, i));
        end % end textInCell
        
        % extract all numbers from cell 'i' in a line of the spreadsheet
        function n = numberInLine(line, i)
            commas = strfind(line, ','); % find indexes of commas
            n = zeros(1, length(commas)-i);
            column = 1;
            for j = i:length(commas)
                str = line(commas(j-1)+1:commas(j)-1);
                if strcmp(str, 'OVER')
                    n(column) = Inf;
                else
                    n(column) = str2double(str);
                end;
                column = column + 1;
            end;
            % get the cells in last columns
            n(column) = str2double(line(commas(j)+1:end));
            % return only array of non NaN
            n = n(~isnan(n));
        end % end textInCell
        
        % filter data using to take out noise
        function dFiltered = filterData(d)
            windowSize = 5;
            dFiltered = filter(ones(1,windowSize)/windowSize,1,d, [], 1);
        end
        
        % calculate specific growth rate of a wavelength
        function r = calculateSpecificRate(time, data)
            % get the number of replicates
            nReps = size(data, 2);
            % log-transform data
            data = log(TimeSeriesData.filterData(data));
            % relicate time array
            timeDifferential = diff(time);
            timeDifferential = repmat(timeDifferential', 1, nReps);
            % filter data to reduce noice
            r = data;
            % calculate derivative
            r = diff(r, 1, 1)./timeDifferential;
            % replicate the first data point to make data same size
            % as original data
            r = [r(1, :); r];
        end;
        
        
        function err = alignmentError(time, curve1, curve2, tau)
            mCurve1 = median(curve1, 2);
            mCurve2 = median(curve2, 2);
            
            indexRef   = find(time<(time(end)-tau));
            indexOther = find(time>tau);
            ref        = mCurve1(indexRef);
            timeRef    = time(indexRef);
            other      = mCurve2(indexOther);
            timeOther  = time(indexOther)-tau;
            
            % take out first time point from timeRef to avoid NaN
            timeRef = timeRef(2:end);
            
            otherInterp = interp1(timeOther, other, timeRef);
            
            ref = ref(2:end);
            otherInterp = otherInterp(:);
            
            
            v = (ref - otherInterp).^2;
            err = sum(v);
        end;
        
        
        function str = num2strRound(n, i)
            
            % str = num2strRound(n, i) - works like num2str but rounds up if next
            % digit is >= 5.
            
            n = n.*10.^i;
            n = round(n);
            n = n./10.^i;
            str = num2str(n);
        end;
        
              % numeric calculation of time derivative. Uses local averaging to
        % smooth data before calculating the finite difference.
        function r = calculateRate(time, data, windowSize)
            % get the number of replicates
            nReps = size(data, 2);
            
            % filter data
            if nargin == 2
                data = TimeSeriesData_jy.filterData2(data);
            else
                data = TimeSeriesData_jy.filterData2(data, windowSize);
            end
            % relicate time array
            timeDifferential = diff(time);
            timeDifferential = repmat(timeDifferential(:), 1, nReps);
            % filter data to reduce noice
            r = data;
            % calculate derivative
            r = diff(r, 1, 1)./timeDifferential;
            % replicate the first data point to make data same size
            % as original data
            r = [r(1, :); r];
        end
        
        % filter data to reduce noise. Optional window size (default: 5)
        function dFiltered = filterData2(d, windowSize)
            if nargin == 1
                windowSize = 5;
            end
            dFiltered = filter(ones(1,windowSize)/windowSize,1,d, [], 1);
        end
        
    end % end static methods
    
    
end % classdef