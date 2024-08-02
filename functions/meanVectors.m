function [meanVector, directVectors] = meanVectors(vectors)
    directVectors = zeros(size(vectors));

    meanVector = vectors(1,:);
    directVectors(1,:) = vectors(1,:);
    for i=2 : size(vectors, 1)
        v1 = vectors(i,:);
        v2 = -v1;

        angle_v1 = acos(dot(v1, meanVector) / (norm(v1) * norm(meanVector)));
        angle_v2 = acos(dot(v2, meanVector) / (norm(v2) * norm(meanVector)));

        if angle_v1 < angle_v2
            v = v1;
        else
            v = v2;
        end

        directVectors(i,:) = v;
        meanVector = meanVector + v;
    end

    meanVector = meanVector ./ size(vectors,1);
end
