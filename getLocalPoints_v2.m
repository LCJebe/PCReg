%% Function: get points within sphere
% returns only the points from pts that are within a certain radius R
% of center point c if there are at least min_points in it
% RETURNS: Points RELATIVE to c
function [pts_sphere, dists] = getLocalPoints_v2(pts, R, c, min_points, max_points)

    % get dists from c
    pts_c = pts-c;
    dists = vecnorm(pts_c, 2, 2);
    
    mask = dists < R;
    num_pts = sum(mask);
    % make sure constraints are met
    if (num_pts < min_points) || (num_pts > max_points)
        dists = [];
        pts_sphere = [];
    else
        dists = dists(find(mask));
        mask = cat(2, mask, mask, mask);
        pts_sphere = reshape(pts_c(mask), [], 3);  
    end
end