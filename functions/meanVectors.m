function [meanVector, directVectors] = meanVectors(vectors)
    % Handle empty input
    if isempty(vectors)
        error('Cannot calculate mean of empty vector array');
    end
    
    % Handle single vector case
    if size(vectors, 1) == 1
        meanVector = vectors(1,:);
        directVectors = vectors;
        return;
    end
    
    directVectors = zeros(size(vectors));

    meanVector = vectors(1,:);
    directVectors(1,:) = vectors(1,:);
    for i=2 : size(vectors, 1)
        v1 = vectors(i,:);
        v2 = -v1;

        % Handle zero vectors to avoid division by zero
        if norm(v1) == 0 || norm(meanVector) == 0
            directVectors(i,:) = v1;
            meanVector = meanVector + v1;
            continue;
        end

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
