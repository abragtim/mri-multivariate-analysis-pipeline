%% Parse inputs
fprintf("Parsing inputs...\n")
settings = readjson("data/settings.json");

patients = {};
patientsFolders = dir('data/');
patientsFolders = patientsFolders(~startsWith({patientsFolders.name}, '.') ...
    & ~strcmp({patientsFolders.name}, 'settings.json') ...
    & ~strcmp({patientsFolders.name}, 'example'));
patientsFoldersCount = length(patientsFolders);
for i = 1:patientsFoldersCount
    patientFolder = patientsFolders(i);
    if(~patientFolder.isdir)
        error('%s is not a patient directory. All objects in the /data directory must be directories.', patientFolder.name);
    end

    intermediateImagesFolder = ['data/', patientFolder.name, '/intermediate_images'];
    if(~exist(intermediateImagesFolder, 'dir'))
        mkdir(intermediateImagesFolder)
    end

    % get metadata
    metadata = readjson("data/" + patientFolder.name + "/metadata.json");
    patientName = patientFolder.name;

    % get T1
    t1PathFolder = ['data/', patientFolder.name, '/'];
    t1 = metric('T1', t1PathFolder, 'T1');
    t1.get();

    % get lesion
    lesionSegmentationPathFolder = ['data/', patientFolder.name, '/'];
    lesionSegmentation = metric('lesionSegmentation', lesionSegmentationPathFolder, 'lesionSegmentation');
    lesionSegmentation.get();

    % get GM
    gmSegmentationPathFolder = ['data/', patientFolder.name, '/'];
    gmSegmentation = metric('GMSegmentation', gmSegmentationPathFolder, 'GMSegmentation');
    gmSegmentation.get();

    % get WM
    wmSegmentationPathFolder = ['data/', patientFolder.name, '/'];
    wmSegmentation = metric('WMSegmentation', wmSegmentationPathFolder, 'WMSegmentation');
    wmSegmentation.get();

    % get metrics
    metrics = cell(length(settings.metrics), 1);
    for j = 1:length(settings.metrics)
        metricName = settings.metrics{j};
        metricPathFolder = ['data/', patientFolder.name, '/'];
        metricEntity = metric(metricName, metricPathFolder, metricName);
        metricEntity.get();
        metrics{j} = metricEntity;
    end

    patients{i} = patient(patientName, metadata.epitype, t1, ...
        lesionSegmentation, wmSegmentation, gmSegmentation, metrics);
    fprintf("Parsed %i/%i patients.\n", i, patientsFoldersCount)
end

%% Transform to MNI. Segment contralateral healthy region in MNI
for i = 1:length(patients)
    patient = patients{i};
    spm_transform_to_mni(patient);  % comment this if MNIs are already computed with SPM

    patient.BrainmaskSegmentationMNI = metric(['mni_', patient.BrainmaskSegmentation.Name], ...
        patient.BrainmaskSegmentation.PathFolder, ...
        ['mni_', patient.BrainmaskSegmentation.NiiFilename]);
    patient.WMSegmentationMNI = metric(['mni_', patient.WMSegmentation.Name], ...
        patient.WMSegmentation.PathFolder, ...
        ['mni_', patient.WMSegmentation.NiiFilename]);
    patient.GMSegmentationMNI = metric(['mni_', patient.GMSegmentation.Name], ...
        patient.GMSegmentation.PathFolder, ...
        ['mni_', patient.GMSegmentation.NiiFilename]);
    patient.LesionSegmentationMNI = metric(['mni_', patient.LesionSegmentation.Name], ...
        patient.LesionSegmentation.PathFolder, ...
        ['mni_', patient.LesionSegmentation.NiiFilename]);

    metricsMNI = cell(length(patient.Metrics), 1);
    for j = 1:length(patient.Metrics)
        patspaceMetric = patient.Metrics{j};
        metricMNIName = ['mni_', patspaceMetric.Name];
        metricMNIFilename = ['mni_', patspaceMetric.NiiFilename];
        metricMNIEntity = metric(metricMNIName, patspaceMetric.PathFolder, metricMNIFilename);
        metricMNIEntity.get();
        metricsMNI{j} = metricMNIEntity;
    end
    patient.MetricsMNI = metricsMNI;

    patient.ContralateralHealthySegmentationMNI = patient.LesionSegmentationMNI.mirror( ...
        'mni_contralateralHealthySegmentation', ...
        [patient.LesionSegmentationMNI.PathFolder, 'intermediate_images/'], ...
        'mni_contralateralHealthySegmentation');

    patients{i} = patient;  % save changes
end

%% Generate results folder
if(~exist('./results', 'dir'))
    mkdir('results')
end

%% Multivariate analysis with CCA
multivariateCoeffsTable = table();
for i = 1:length(patients)
    patient = patients{i};
    fprintf("Performing multivariate analysis for patient %s...\n", patient.Name)

    brainmaskSegment = logical(patient.BrainmaskSegmentationMNI.get().img);
    lesionSegment = logical(patient.LesionSegmentationMNI.get().img) & brainmaskSegment;
    contralateralSegment = logical(patient.ContralateralHealthySegmentationMNI.get().img) & brainmaskSegment;

    dataForCCA = [];
    for j = 1:length(patient.MetricsMNI)
        metric = patient.MetricsMNI{j};
        metricImg = metric.get().img;
        lesioned = metricImg(lesionSegment);
        healthy = metricImg(contralateralSegment);
        dataForCCA = [dataForCCA, [lesioned; healthy]];
    end
    labelsForCCA = logical([ones(size(lesioned)); zeros(size(healthy))]);

    dataNormalizedCCA = zscore(dataForCCA, 0, 1);
    [coeffs, ~, correlation_r, ~, ~] = canoncorr(dataNormalizedCCA, labelsForCCA);

    patient.ContralateralCCAResult = ccaResult(coeffs, correlation_r);
    patients{i} = patient;  % save changes


    % visualize results
    resultsFolder = [metric.PathFolder, '/results'];
    if(~exist(resultsFolder, 'dir'))
        mkdir(resultsFolder)
    end

    metricNames = cellfun(@(m) ['coeff_', m.Name], patient.MetricsMNI, 'UniformOutput', false);
    newRow = cell2table([{string(patient.Name)}, num2cell(coeffs)', num2cell(correlation_r)'], ...
                   'VariableNames', [{'Patient'}, metricNames', {'CorrelationR'}]);
    multivariateCoeffsTable = [multivariateCoeffsTable; newRow];

    fig = figure('Visible','on');
    spider_plot(coeffs', 'AxesLimits',...
        [min(coeffs) .* ones(size(coeffs))'; max(coeffs) .* ones(size(coeffs))'], ...
        'AxesLabels', cellfun(@(name) strrep(name, '_', ' '), metricNames,...
        'UniformOutput', false))
    title(['Contralateral CCA ', patient.Name])
    savefig(fig, [resultsFolder, '/spider_plot_cca_coeffs.fig'])
end
disp(multivariateCoeffsTable)
save('results/multivariateCoeffsTable.mat', 'multivariateCoeffsTable');

%% Calculate and generate new CCA metrics
ccaMeanResult = calculateCCAmeanResult(patients, settings.correlation_r_significant_threshold);

% mean CCA per epitype
patientsEpitypes = cellfun(@(p) p.Epitype, patients, 'UniformOutput', false);
uniqueEpitypes = unique(patientsEpitypes);
ccaEpitypesMeanResults = cell(size(uniqueEpitypes));
for i = 1:length(uniqueEpitypes)
    epitype = uniqueEpitypes{i};
    patientsWithEpitype = strcmp(patientsEpitypes, epitype);
    ccaEpitypesMeanResults{i} = calculateCCAmeanResult(patients(patientsWithEpitype), ...
        settings.correlation_r_significant_threshold);
end
ccaMeanResultPerEpitype = containers.Map(uniqueEpitypes, ccaEpitypesMeanResults);

% generate new metrics
for i = 1:length(patients)
    patient = patients{i};

    % combined metric with maximum correlation r for this patient
    ccaMaxCorrMetric = patient.ContralateralCCAResult.calculateCombinedMetric( ...
        patient.MetricsMNI, patient.BrainmaskSegmentationMNI, 'ccaMaxCorr');
    ccaMaxCorrMetric.equalize(patient.BrainmaskSegmentationMNI);

    ccaMeanMetric = ccaMeanResult.calculateCombinedMetric(patient.MetricsMNI, patient.BrainmaskSegmentationMNI, 'ccaMean');
    ccaMeanMetric.equalize(patient.BrainmaskSegmentationMNI);

    ccaMeanPerEpitypeMetric = ccaMeanResultPerEpitype(patient.Epitype)...
        .calculateCombinedMetric(patient.MetricsMNI, patient.BrainmaskSegmentationMNI, ['ccaMean', '_', patient.Epitype]);
    ccaMeanPerEpitypeMetric.equalize(patient.BrainmaskSegmentationMNI)

    patient.MetricsMNI = [patient.MetricsMNI(:)', ...
        {ccaMaxCorrMetric}, {ccaMeanMetric}, {ccaMeanPerEpitypeMetric}];
    patients{i} = patient;  % save changes
end

%% Univariate analysis
univariateAnalysisTable = table();
for i = 1:length(patients)
    patient = patients{i};

    brainmaskSegmentation = logical(patient.BrainmaskSegmentationMNI.get().img);
    wmSegmentation = logical(patient.WMSegmentationMNI.get().img);
    gmSegmentation = logical(patient.GMSegmentationMNI.get().img);
    lesionSegmentation = logical(patient.LesionSegmentationMNI.get().img);
    contralateralSegmentation = logical(patient.ContralateralHealthySegmentationMNI.get().img);

    for j = 1:length(patient.MetricsMNI)
        metric = patient.MetricsMNI{j};
        [meanLesion, stdLesion, medianLesion, ...
            meanHealthy, stdHealthy, medianHealthy, ...
            meanHealthyWM, stdHealthyWM, medianHealthyWM, ...
            meanHealthyGM, stdHealthyGM, medianHealthyGM, ...
            meanDiff, stdDiff, medianDiff, ...
            pval, zval, correlationR] = calculateUnivariateStatistics( ...
            metric, lesionSegmentation & brainmaskSegmentation, ...
            contralateralSegmentation & brainmaskSegmentation);

        newRow = table(string(patient.Name), string(patient.Epitype), ...
            string(metric.Name), meanLesion, stdLesion, medianLesion, ...
            meanHealthy, stdHealthy, medianHealthy, ...
            meanHealthyWM, stdHealthyWM, medianHealthyWM, ...
            meanHealthyGM, stdHealthyGM, medianHealthyGM, ...
            meanDiff, stdDiff, medianDiff,...
            pval, zval, correlationR, ...
            'VariableNames', {'Patient', 'Epitype', 'Metric', 'MeanLesion', 'StdLesion', 'MedianLesion', ...
            'MeanHealthy', 'StdHealthy', 'MedianHealthy', ...
            'MeanHealthyWM', 'StdHealthyWM', 'MedianHealthyWM', ...
            'MeanHealthyGM', 'StdHealthyGM', 'MedianHealthyGM', ...
            'MeanPercentageDiff', 'StdPercentageDiff', 'MedianPercentageDiff',...
            'PValue', 'ZValue', 'CorrelationR'});
        univariateAnalysisTable = [univariateAnalysisTable; newRow];
    end
end
disp(univariateAnalysisTable)
save('results/univariateAnalysisTable.mat', 'univariateAnalysisTable');

% visualize
function [] = saveBoxchart(metrics, values, epitypes, label)
    fig = figure('Visible','on');
    epitypes = categorical(strrep(epitypes, '_', ' '));
    boxchart(categorical(strrep(metrics, '_', ' ')), values, ...
        'GroupByColor', epitypes);
    ylabel(label)
    legend(unique(epitypes))
    savefig(fig, ['results/boxchart_', lower(strrep(label, ' ', '_')), '.fig'])
end

saveBoxchart(univariateAnalysisTable.Metric,...
    univariateAnalysisTable.PValue,...
    univariateAnalysisTable.Epitype, ...
    'p-value');
saveBoxchart(univariateAnalysisTable.Metric,...
    univariateAnalysisTable.CorrelationR,...
    univariateAnalysisTable.Epitype, ...
    'Correlation R');
saveBoxchart(univariateAnalysisTable.Metric,...
    univariateAnalysisTable.MeanPercentageDiff,...
    univariateAnalysisTable.Epitype, ...
    'Mean percentage difference');
saveBoxchart(univariateAnalysisTable.Metric,...
    univariateAnalysisTable.StdPercentageDiff,...
    univariateAnalysisTable.Epitype, ...
    'SD percentage difference');
saveBoxchart(univariateAnalysisTable.Metric,...
    univariateAnalysisTable.MedianPercentageDiff,...
    univariateAnalysisTable.Epitype, ...
    'Median percentage difference');

fprintf('Completed successfully!\n')
