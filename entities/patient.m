classdef patient
    properties
        Name
        Epitype
        T1 metric

        LesionSegmentation metric
        LesionSegmentationMNI metric

        ContralateralHealthySegmentationMNI metric

        BrainmaskSegmentation metric
        BrainmaskSegmentationMNI metric

        WMSegmentation
        WMSegmentationMNI

        GMSegmentation
        GMSegmentationMNI

        Metrics cell
        MetricsMNI cell

        ContralateralCCAResult ccaResult
    end

    methods
        function obj = patient(name, epitype, t1, lesionSegmentation, ...
                WMSegmentation, GMSegmentation, metrics)
            obj.Name = name;
            obj.Epitype = epitype;
            obj.T1 = t1;
            obj.LesionSegmentation = lesionSegmentation;
            obj.Metrics = metrics;
            obj.WMSegmentation = WMSegmentation;
            obj.GMSegmentation = GMSegmentation;

            obj.BrainmaskSegmentation = obj.calculateBrainmask();
        end
    end

    methods(Access=private)
        function brainmask = calculateBrainmask(obj)
            wmNii = obj.WMSegmentation.get();
            gmNii = obj.GMSegmentation.get();

            brainmaskNii = wmNii;
            brainmaskNii.img = wmNii.img | gmNii.img;  % todo: test it

            brainmask = metric('brainmask', ...
                [obj.WMSegmentation.PathFolder, 'intermediate_images/'], ...
                'brainmask');
            brainmask.save(brainmaskNii);
        end    
    end
end
