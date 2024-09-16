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

    hold off;
    savefig(fig, [metric.PathFolder, '/results', '/hist_', metric.Name, '.fig'])
    close(fig);

    %% single-value metrics
    meanLesion = mean(lesion);
    medianLesion = median(lesion);
    stdLesion = std(lesion);

    meanContralateralHealthy = mean(healthy);
    medianContralateralHealthy = median(healthy);
    stdContralateralHealthy = std(healthy);

    meanContralateralHealthyWM = mean(healthyWM);
    medianContralateralHealthyWM = median(healthyWM);
    stdContralateralHealthyWM = std(healthyWM);

    meanContralateralHealthyGM = mean(healthyGM);
    medianContralateralHealthyGM = median(healthyGM);
    stdContralateralHealthyGM = std(healthyGM);

    meanPercentageDiff = (meanLesion - meanContralateralHealthy) / meanContralateralHealthy * 100;
    medianPercentageDiff = (medianLesion - medianContralateralHealthy) / medianContralateralHealthy * 100;
    stdPercentageDiff = (stdLesion - stdContralateralHealthy) / stdContralateralHealthy * 100;
end
