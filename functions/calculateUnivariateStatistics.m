function [meanLesion, stdLesion, medianLesion, ...
    meanContralateralHealthy, stdContralateralHealthy, medianContralateralHealthy,...
    meanPercentageDiff, stdPercentageDiff, medianPercentageDiff,...
    pval, zval, correlationR] ...
    = calculateUnivariateStatistics(metric, lesionSegmentation, contralateralSegmentation)

    volume = metric.get().img;
    lesion = volume(lesionSegmentation);
    healthy = volume(contralateralSegmentation);

    %% U - test
    [pval, ~, stats] = ranksum(lesion, healthy);
    zval = stats.zval;
    correlationR = abs(zval) / sqrt(numel(lesion) + numel(healthy));

    %% save hist
    fig = figure('Visible','on');
    hold on;
    histPositive = histogram(lesion);
    lesionCounts = histPositive.NumBins;
    healthyCounts = floor((max(healthy) - min(healthy)) ./ (max(lesion) - min(lesion)) .* lesionCounts);
    histogram(healthy, 'NumBins', healthyCounts);
    ylabel('Counts (-)')
    legend('Lesion', 'Contralateral healthy region')

    resultsFolder = [metric.PathFolder, '/results'];
    if(~exist(resultsFolder, 'dir'))
        mkdir(resultsFolder)
    end
    hold off;
    savefig(fig, [resultsFolder, '/hist_', metric.Name, '.fig'])
    close(fig);

    %% single-value metrics
    meanLesion = mean(lesion);
    medianLesion = median(lesion);
    stdLesion = std(lesion);

    meanContralateralHealthy = mean(healthy);
    medianContralateralHealthy = median(healthy);
    stdContralateralHealthy = std(healthy);

    meanPercentageDiff = (meanLesion - meanContralateralHealthy) / meanContralateralHealthy * 100;
    medianPercentageDiff = (medianLesion - medianContralateralHealthy) / medianContralateralHealthy * 100;
    stdPercentageDiff = (stdLesion - stdContralateralHealthy) / stdContralateralHealthy * 100;
end
