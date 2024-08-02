function [ccaMeanResult] = calculateCCAmeanResult(patients, correlation_r_significant_threshold)
    ccaResults = cellfun(@(p) p.ContralateralCCAResult, patients, 'UniformOutput', false);
    coeffs = cellfun(@(x) x.Coefficients, ccaResults, 'UniformOutput', false);
    coeffs = cell2mat(coeffs)';
    correlationRs = cellfun(@(x) x.MaxCorrelationR, ccaResults, 'UniformOutput', false);
    correlationRs = cell2mat(correlationRs)';

    withEffect = correlationRs > correlation_r_significant_threshold;
    coeffsWithEffect = coeffs(withEffect, :);

    [ccaMeanCoeffs, ~] = meanVectors(coeffsWithEffect);
    meanCorrelationR = mean(correlationRs(withEffect));
    ccaMeanResult = ccaResult(ccaMeanCoeffs, meanCorrelationR);
end
