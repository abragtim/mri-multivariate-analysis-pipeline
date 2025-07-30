function [meanLesion, stdLesion, medianLesion, ...
    meanContralateralHealthy, stdContralateralHealthy, medianContralateralHealthy,...
    meanContralateralHealthyWM, stdContralateralHealthyWM, medianContralateralHealthyWM,...
    meanContralateralHealthyGM, stdContralateralHealthyGM, medianContralateralHealthyGM,...
    meanPercentageDiff, stdPercentageDiff, medianPercentageDiff,...
    pval, zval, correlationR] ...
    = calculateUnivariateStatistics(metric, lesionSegmentation, contralateralSegmentation, ...
                                    brainmask, wmSegmentation, gmSegmentation)

    volume = metric.get().img;
    lesion = volume(lesionSegmentation & brainmask);
    healthy = volume(contralateralSegmentation & brainmask);
    healthyWM = volume(contralateralSegmentation & wmSegmentation);
    healthyGM = volume(contralateralSegmentation & gmSegmentation);

    %% U - test
    % Check if we have sufficient data for statistical test
    if isempty(lesion) || isempty(healthy) || length(lesion) < 2 || length(healthy) < 2
        fprintf('Warning: Insufficient data for statistical test in metric %s. Using default values.\n', metric.Name);
        pval = 1.0;  % No significant difference
        zval = 0;
        correlationR = 0;
    else
        try
            [pval, ~, stats] = ranksum(lesion, healthy);
            zval = stats.zval;
            correlationR = abs(zval) / sqrt(numel(lesion) + numel(healthy));
        catch ME
            fprintf('Warning: Statistical test failed for metric %s: %s. Using default values.\n', metric.Name, ME.message);
            pval = 1.0;
            zval = 0;
            correlationR = 0;
        end
    end

    %% save hist
    try
        fig = figure('Visible','on');
        hold on;
        
        % Check if we have valid data for histogram
        if ~isempty(lesion) && length(unique(lesion)) > 1
            histPositive = histogram(lesion);
            lesionCounts = histPositive.NumBins;
        else
            lesionCounts = 10;  % Default bin count
        end
        
        if ~isempty(healthy) && length(unique(healthy)) > 1
            % Calculate healthy bin count safely
            if ~isempty(lesion) && max(lesion) > min(lesion) && max(healthy) > min(healthy)
                healthyCounts = max(1, floor((max(healthy) - min(healthy)) ./ (max(lesion) - min(lesion)) .* lesionCounts));
            else
                healthyCounts = lesionCounts;
            end
            histogram(healthy, 'NumBins', healthyCounts);
        end
        
        ylabel('Counts (-)')
        legend('Lesion', 'Contralateral healthy region')
        hold off;
        
        savefig(fig, [metric.PathFolder, '/results', '/hist_', metric.Name, '.fig'])
        close(fig);
    catch ME
        fprintf('Warning: Histogram visualization failed for metric %s: %s\n', metric.Name, ME.message);
        if exist('fig', 'var') && isvalid(fig)
            close(fig);
        end
    end

    %% single-value metrics
    % Handle empty arrays gracefully
    if isempty(lesion)
        meanLesion = NaN;
        medianLesion = NaN;
        stdLesion = NaN;
    else
        meanLesion = mean(lesion);
        medianLesion = median(lesion);
        stdLesion = std(lesion);
    end

    if isempty(healthy)
        meanContralateralHealthy = NaN;
        medianContralateralHealthy = NaN;
        stdContralateralHealthy = NaN;
    else
        meanContralateralHealthy = mean(healthy);
        medianContralateralHealthy = median(healthy);
        stdContralateralHealthy = std(healthy);
    end

    if isempty(healthyWM)
        meanContralateralHealthyWM = NaN;
        medianContralateralHealthyWM = NaN;
        stdContralateralHealthyWM = NaN;
    else
        meanContralateralHealthyWM = mean(healthyWM);
        medianContralateralHealthyWM = median(healthyWM);
        stdContralateralHealthyWM = std(healthyWM);
    end

    if isempty(healthyGM)
        meanContralateralHealthyGM = NaN;
        medianContralateralHealthyGM = NaN;
        stdContralateralHealthyGM = NaN;
    else
        meanContralateralHealthyGM = mean(healthyGM);
        medianContralateralHealthyGM = median(healthyGM);
        stdContralateralHealthyGM = std(healthyGM);
    end

    % Handle division by zero in percentage calculations
    if isnan(meanContralateralHealthy) || meanContralateralHealthy == 0
        meanPercentageDiff = NaN;
    else
        meanPercentageDiff = (meanLesion - meanContralateralHealthy) / meanContralateralHealthy * 100;
    end
    
    if isnan(medianContralateralHealthy) || medianContralateralHealthy == 0
        medianPercentageDiff = NaN;
    else
        medianPercentageDiff = (medianLesion - medianContralateralHealthy) / medianContralateralHealthy * 100;
    end
    
    if isnan(stdContralateralHealthy) || stdContralateralHealthy == 0
        stdPercentageDiff = NaN;
    else
        stdPercentageDiff = (stdLesion - stdContralateralHealthy) / stdContralateralHealthy * 100;
    end
end
