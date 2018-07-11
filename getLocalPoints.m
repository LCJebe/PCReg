%% Function: get points within sphere
% returns only the points from pts that are within a certain radius R
% of center point c if there are at least min_points in it
% RETURNS: Points RELATIVE to c, optinally sorted by distance. 
function [pts_sphere, dists] = getLocalPoints(pts, R, c, min_points)

    % first: quickly crop cube around the center
    xLim = c(1) + [-R, R];
    yLim = c(2) + [-R, R];
    zLim = c(3) + [-R, R];
    mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
         & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
         & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
    mask = cat(2, mask, mask, mask);
    pts_cube = reshape(pts(mask), [], 3);    
    
    if size(pts_cube, 1) < min_points
        pts_sphere = [];
        dists = [];
    else
        % then: get distances for each remaining point and only keep the points
        % that are within the radius
        pts_rel = pts_cube - c;
        dists = vecnorm(pts_rel, 2, 2);
        mask = dists < R;
        dists = dists(find(mask));
        mask = cat(2, mask, mask, mask);
        pts_sphere = reshape(pts_rel(mask), [], 3);    
        
        % reject if not enough points found
        if size(pts_sphere, 1) < min_points
            pts_sphere = [];
            dists = [];
        end
    end
end