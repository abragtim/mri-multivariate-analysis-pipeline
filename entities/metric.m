classdef metric
    properties
        Name
        PathFolder
        NiiFilename
    end

    methods
        function obj = metric(name, pathFolder, niiFilename)
            obj.Name = name;
            obj.PathFolder = pathFolder;
            obj.NiiFilename = niiFilename;
        end

        function niiImage = get(obj)
            niiImage = load_untouch_nii(obj.getPath());
        end

        function save(obj, nii)
            save_untouch_nii(nii, obj.getPath())
        end

        function mirrored = mirror(obj, mirroredName, mirroredPathFolder, mirroredNiiFilename)
            originalNii = obj.get();
            mirroredNii = originalNii;
            mirroredNii.img = flip(originalNii.img, 1);

            mirrored = metric(mirroredName, mirroredPathFolder, mirroredNiiFilename);
            mirrored.save(mirroredNii);
        end

        function equalized = equalize(obj, brainmask)
            brainmaskSegment = logical(brainmask.get().img);
            originalNii = obj.get();
            equalizedNii = originalNii;
            equalizedNii.img(brainmaskSegment) = histeq(originalNii.img(brainmaskSegment));
            equalizedNii.img = applyWindowing(equalizedNii.img, brainmaskSegment, 1, 99);

            equalized = metric([obj.Name, '_equalized'], obj.PathFolder, [obj.NiiFilename, '_equalized']);
            equalized.save(equalizedNii);
        end

        function path = getPath(obj)
            path = [obj.PathFolder, obj.NiiFilename, '.nii'];
        end
    end
end
