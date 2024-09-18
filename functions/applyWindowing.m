function [res] = applyWindowing(img,brainmask, low, high)
    p_low = prctile(img(brainmask), low);
    p_high = prctile(img(brainmask), high);
    res = img;
    res(brainmaskSegment & res < p_low) = p_low;
    res(brainmaskSegment & res > p_high) = p_high;
    res(~brainmaskSegment) = min(res(brainmask), [], 'all');
end
