%% function to calculate a descriptor for each point
function [feat, desc] = getRotationalDescriptors(pts, sample_pts, options)
    % pts: points in pointcloud
    % sample_pts: points to calculate descriptors at
    % min_points: minimum number of points in sphere
    % R: Radius of sphere
    % thVar: two element vector that contains the two thresholds for the
    % eigenvalues of the covariance matrix (sphere-reject)
    
    % returns:
        % - feat: feature locations
        % - desc: feature descriptors
        
    min_pts = options.min_pts;
    R = options.R;

    % determine phases to be used for standing wave descriptor
    phases = 0:pi/2:pi*4;    
            
    % preallocate space for features and descriptors, even though we don't
    % know their length    
    desc = nan(size(sample_pts, 1), 2*length(phases)-1);
    feat = nan(size(sample_pts, 1), 3);
    tic
    for i = 1:size(sample_pts, 1)
        c = sample_pts(i, :);
        
        % return local points (only care about distances)
        [~, dists] = getLocalPoints(pts, R, c, min_pts);

        if ~ isempty(dists) 
            
            % ---- calculate the descriptor
            desc_entry = zeros(1, 2*length(phases)-1);
            desc_entry(1) = sum(1*dists);
            for j = 2:length(phases)
               ph = phases(j);
               desc_entry(j) = sum(sin(ph*dists));
               desc_entry(j+1) = sum(cos(ph*dists));
            end

            % descriptor structure is 
            % [r, sin(r*0.5pi), cos(r*0.5pi), sin(r*pi), cos(r*pi),
            % sin(r*1.5pi), cos(r*15.pi), ...]

            desc(i, :) = desc_entry;
            feat(i, :) = c;
        end
    end
    
    % remove nan rows from desc and feat and return
    mask = find(~isnan(desc(:, 1))); % row indices of filled rows
    desc = desc(mask, :);
    feat = feat(mask, :);
    toc
end
