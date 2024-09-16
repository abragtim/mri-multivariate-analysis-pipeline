function [outputArg1,outputArg2] = visualizeAsBoxchart(inputArg1,inputArg2)
    fig = figure('Visible','on');
    boxchart(categorical(strrep(univariateAnalysisTable.Metric, '_', ' ')),...
        univariateAnalysisTable.CorrelationR, ...
        'GroupByColor', categorical(strrep(univariateAnalysisTable.Epitype, '_', ' ')));
    ylabel('Correlation R')
    if(~exist('./results', 'dir'))
        mkdir('results')
    end
    savefig(fig, 'results/boxchart_correlation_r.fig')
end

