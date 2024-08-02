function spm_transform_to_mni(patientEntity)
arguments
    patientEntity patient
end
    jobfile = {'functions/spm_transform_to_mni_job.m'};
    jobs = repmat(jobfile, 1, 1);
    inputs = cell(2, 1);
    inputs{1, 1} = {patientEntity.T1.getPath()}; % Normalise: Estimate & Write: Image to Align - cfg_files

    metricsCounts = length(patientEntity.Metrics);
    inputsMetrics = cell(metricsCounts + 2, 1); % metrics + brainmask + lesion
    for i = 1:metricsCounts
        inputsMetrics{i} = patientEntity.Metrics{i}.getPath();
    end
    inputsMetrics{metricsCounts + 1} = patientEntity.BrainmaskSegmentation.getPath();
    inputsMetrics{metricsCounts + 2} = patientEntity.LesionSegmentation.getPath();
    inputs{2, 1} = inputsMetrics;

    spm('defaults', 'PET');
    spm_jobman('run', jobs, inputs{:});
end
