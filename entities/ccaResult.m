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
            
            brainmaskedNormalized = cell(size(metrics));
            for i=1:length(metrics)
                metricImg = metrics{i}.get().img;
                brainmaskedNormalized{i} = metricImg;
                brainmaskedNormalized{i}(brainmaskSegment) = zscore(metricImg(brainmaskSegment));
            end

            metricsSegments = cat(4, brainmaskedNormalized{:});
            combinedMetricImg = sum(bsxfun(@times, metricsSegments, reshape(obj.Coefficients, 1, 1, 1, length(obj.Coefficients))), 4);
            combinedMetricImg(~brainmaskSegment) = min( ...
                combinedMetricImg(brainmaskSegment), [], 'all');
            combinedMetricImg = mat2gray(combinedMetricImg);

            combinedMetricNii = metrics{1}.get();
            combinedMetricNii.img = combinedMetricImg;
            combinedMetric = metric(metricName, metrics{1}.PathFolder, metricName);
            combinedMetric.save(combinedMetricNii);
        end
    end
end
