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
            combinedMetricImg = combinedMetricImg ./ norm(obj.Coefficients);

            % windowing
            p99 = prctile(combinedMetricImg(brainmaskSegment), 99);
            p1 = prctile(combinedMetricImg(brainmaskSegment), 1);
            combinedMetricImg(brainmaskSegment & combinedMetricImg < p1) = p1;
            combinedMetricImg(brainmaskSegment & combinedMetricImg > p99) = p99;
            combinedMetricImg(~brainmaskSegment) = min(combinedMetricImg(brainmaskSegment), [], 'all');

            combinedMetricNii = metrics{1}.get();
            combinedMetricNii.img = combinedMetricImg;
            combinedMetric = metric(metricName, metrics{1}.PathFolder, metricName);
            combinedMetric.save(combinedMetricNii);
        end
    end
end
