function spm_transform_to_patspace(patientEntity)
arguments
    patientEntity patient
end
    inputs = cell(2, 1);
    inputs{1, 1} = {[patientEntity.T1.PathFolder, '/y_T1.nii']};
    inputs{2, 1} = {[patientEntity.T1.PathFolder, '/T1.nii']};

    originalCounts = length(patientEntity.Metrics);
    metricsCounts = length(patientEntity.MetricsMNI) - originalCounts;
    inputsMetrics = cell(metricsCounts * 2, 1);  % CCA metrics + equalized
    for i = 1:metricsCounts
        metric = patientEntity.MetricsMNI{originalCounts + i};
        inputsMetrics{i} = metric.getPath();
        inputsMetrics{metricsCounts + i} = [metric.PathFolder, metric.NiiFilename, '_equalized.nii'];
    end
    inputs{3, 1} = inputsMetrics;

    spm('defaults', 'PET');
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = inputs{1,1};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = inputs{2,1};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = inputs{3,1};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = {patientEntity.T1.PathFolder};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'patspace_';
    spm_jobman('run', matlabbatch);
end
