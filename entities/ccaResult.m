classdef ccaResult    
    properties
        Coefficients
        MaxCorrelationR
    end
    
    methods
        function obj = ccaResult(coeffs, r)
            obj.Coefficients = coeffs;
            obj.MaxCorrelationR = r;
        end

        function combinedMetric = calculateCombinedMetric(obj, metrics, brainmask, metricName)
            brainmaskSegment = logical(brainmask.get().img);
            
            metricsSegments = cellfun(@(m) m.get().img, ...
                metrics, 'UniformOutput',false);
            metricsSegments = cat(4, metricsSegments{:});
            metricsSegments(brainmaskSegment) = zscore(metricsSegments(brainmaskSegment), 0, 4);
            combinedMetricImg = sum(bsxfun(@times, metricsSegments, reshape(obj.Coefficients, 1, 1, 1, length(obj.Coefficients))), 4);

            combinedMetricNii = metrics{1}.get();
            combinedMetricNii.img = combinedMetricImg;
            combinedMetric = metric(metricName, metrics{1}.PathFolder, metricName);
            combinedMetric.save(combinedMetricNii);
        end
    end
end
