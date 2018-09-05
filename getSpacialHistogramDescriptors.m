%% function to calculate a descriptor for each point
function [feat, desc] = getSpacialHistogramDescriptors(pts, sample_pts, options)
    % pts: points in pointcloud
    % sample_pts: points to calculate descriptors at
    % options.min_pts: minimum number of points in sphere
    % options.R: Radius of sphere
    % options.thVar: two element vector that contains the two thresholds for the
        % eigenvalues of the covariance matrix (sphere-reject)
    % options.ALIGN_POINTS: use local reference frame?
    % options.CENTER: center to centroid before calculating descriptor?
    % options.k: k nearest neighbors used for alignment
    
    % returns:
        % - feat: feature locations
        % - desc: feature descriptors
        
    % unpack options
    min_pts = options.min_pts;
    max_pts = options.max_pts;
    R = options.R;
    thVar = options.thVar;
    K = options.k;
    ALIGN_POINTS = options.ALIGN_POINTS;
    CENTER = options.CENTER;
    if isfield(options, 'VERBOSE')
        VERBOSE = options.VERBOSE;
    else
        VERBOSE = 1;
    end
        
    % more internal options
    NORMALIZE = false; % should be true unless point density is similar    
    ALIGN_v2 = true; % use new method for alignment to see if it performs better. 
    r = 2.5; % radius for align_v2.
        
    % define histogram layout (bins) in spherical coordinates
    % it should be NUM_PHI = 2*NUM_THETA
    NUM_R = 10; % 10
    NUM_THETA = 7; % 7
    NUM_PHI = 14; % 14
    
    % in a first step, throw out keypoints that don't fullfill the
    % min_points and max_points constraints
    if VERBOSE
        tic
    end
    valid_mask = false(size(sample_pts, 1), 1);
    parfor i = 1:size(sample_pts, 1)
        c = sample_pts(i, :);
        [~, dists] = getLocalPoints(pts, R, c, min_pts, max_pts);
        if ~isempty(dists)
            valid_mask(i) = 1;
        end
    end
    
    % retain only valid sample_pts
    valid_mask = cat(2, valid_mask, valid_mask, valid_mask);
    sample_pts = reshape(sample_pts(valid_mask), [], 3);
       
    % preallocate space for features and descriptors    
    desc = nan(size(sample_pts, 1), NUM_R * NUM_PHI * NUM_THETA);
    feat = nan(size(sample_pts, 1), 3);
    
    parfor i = 1:size(sample_pts, 1) % PARFOR
        c = sample_pts(i, :);
        
        % return local points
        [pts_local, ~] = getLocalPoints(pts, R, c, min_pts, max_pts);
        
        if ~ isempty(pts_local) 
            num_points = size(pts_local, 1);
            
            % ---- rejection based on PCA variances (== eigenvalues of
            % covariance matrix!!) of k nearest neighbors  
            if strcmp(K, 'all') || K == 1
                k = num_points;
            else
                k = round(num_points*K);
                % sort points by distance to centroid for KNN
                centroid = mean(pts_local, 1);
                dists = vecnorm(pts_local - centroid, 2, 2);
                [~, I] = sort(dists);
                pts_local = pts_local(I, :);
            end
            pts_k = pts_local(1:k, :);
            if ~(sum(thVar == 1) == 2) || ALIGN_POINTS
                % use old alignment method
                if ~ALIGN_v2
                    %[coeff, pts_lrf, variances] = pca(pts_k, 'Algorithm', 'eig', 'Centered', false);
                    
                    [coeff, pts_lrf, variances] = pca(pts_k, 'Algorithm', 'eig');
                % use new method for alignment
                else 
                    % --- 1) find median (L1) / mean (L2)
                    centroid = mean(pts_k, 1);

                    % --- 2) get local neighborhood around this point with radius r and
                    % make sure there is at least a bare minimum of points included
                    [pts_rel, ~] = getLocalPoints(pts_k, r, centroid, 250, inf);

                    if isempty(pts_rel)
                        continue
                    else
                        % --- 3) use those points for alignment (get transform)
                        %[coeff, pts_lrf, variances] = pca(pts_rel, 'Algorithm', 'eig', 'Centered', false); 
                        [coeff, pts_lrf, variances] = pca(pts_rel, 'Algorithm', 'eig'); 
                    end
                end
                
                % continue if constraints are not met
                if (variances(1) / variances(2) < thVar(1)) || ...
                        (variances(2) / variances(3) < thVar(2))
                    continue
                end
            end

            % ---- use sign disambiguition method for aligned points
            k = size(pts_lrf, 1);
            % count number of points with positive sign and see if they
            % dominate ( k nearest neighbors)
            if ALIGN_POINTS
                x_sign = sum(sign(pts_lrf(:, 1)) == 1) >= k/2;
                z_sign = sum(sign(pts_lrf(:, 3)) == 1) >= k/2;

                % map from {0, 1} to {-1, 1}
                x_sign = x_sign*2-1;
                z_sign = z_sign*2-1;

                %  get y sign so that rotation is proper
                y_sign = det(coeff .* [x_sign, 1, z_sign]);
                
                % apply signs to transform matrix (pca coefficients)
                coeff_unambig = coeff .* [x_sign, y_sign, z_sign];
                
                % transform all points into the new coordinate system
                if CENTER
                    pts_local = (pts_local - mean(pts_local, 1))*coeff_unambig;
                else
                    pts_local = pts_local*coeff_unambig;
                end
            end
            
            % ---- calculate the descriptor
           
            % transform points into spherical coordinates
            r_spheric = vecnorm(pts_local, 2, 2);
            theta_spheric = acos(pts_local(:, 3) ./ r_spheric);
            phi_spheric = atan2(pts_local(:, 2), pts_local(:, 2));
            pts_spheric = [r_spheric, theta_spheric, phi_spheric];
            
            r_equi = 0:R^3/NUM_R:R^3;
            r_bins = nthroot(r_equi, 3);
            phi_bins = -pi:2*pi/NUM_PHI:pi;
            theta_bins = 0:pi/NUM_THETA:pi;
            
            % get trivariate histogram
            [counts, ~, ~, ~] = histcn(pts_spheric, r_bins, theta_bins, phi_bins);
            
            % flatten 3D histogram to 1D
            new_entry = reshape(counts, [], 1);
            
            % optional: normalize
            if NORMALIZE
                new_entry = new_entry / size(pts_local, 1);
            end
            
            desc(i, :) = new_entry;
            feat(i, :) = c;
        end
    end
    
    % remove nan rows from desc and feat and return
    mask = find(~isnan(desc(:, 1))); % row indices of filled rows
    desc = desc(mask, :);
    feat = feat(mask, :);
    if VERBOSE
        fprintf('Calculated descriptors in %0.1f seconds...\n', toc);
    end
end
