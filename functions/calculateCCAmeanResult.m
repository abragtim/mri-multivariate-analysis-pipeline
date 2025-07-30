function [ccaMeanResult] = calculateCCAmeanResult(patients, correlation_r_significant_threshold)
    % Handle empty patients array
    if isempty(patients)
        error('No patients provided for CCA mean calculation');
    end
    
    ccaResults = cellfun(@(p) p.ContralateralCCAResult, patients, 'UniformOutput', false);
    coeffs = cellfun(@(x) x.Coefficients, ccaResults, 'UniformOutput', false);
    
    % Handle single patient case - ensure proper matrix dimensions
    if length(coeffs) == 1
        coeffs = coeffs{1}(:)';  % Ensure row vector
    else
        coeffs = cell2mat(coeffs)';
    end
    
    correlationRs = cellfun(@(x) x.MaxCorrelationR, ccaResults, 'UniformOutput', false);
    
    % Handle single patient case for correlation Rs
    if length(correlationRs) == 1
        correlationRs = correlationRs{1};
    else
        correlationRs = cell2mat(correlationRs)';
    end

    withEffect = correlationRs > correlation_r_significant_threshold;
    
    % Handle case where no patients meet the threshold
    if ~any(withEffect)
        fprintf('Warning: No patients meet the correlation R threshold (%.3f). Using all patients.\n', correlation_r_significant_threshold);
        withEffect = true(size(correlationRs));
    end
    
    coeffsWithEffect = coeffs(withEffect, :);

    [ccaMeanCoeffs, ~] = meanVectors(coeffsWithEffect);
    meanCorrelationR = mean(correlationRs(withEffect));
    ccaMeanResult = ccaResult(ccaMeanCoeffs, meanCorrelationR);
end
