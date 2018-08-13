%% function to calculate a descriptor for each point
function [feat, desc] = getSpacialHistogramDescriptors(pts, sample_pts, min_pts, R, thVar, ALIGN_POINTS, CENTER, K)
    % pts: points in pointcloud
    % sample_pts: points to calculate descriptors at
    % min_points: minimum number of points in sphere
    % R: Radius of sphere
    % thVar: two element vector that contains the two thresholds for the
    % eigenvalues of the covariance matrix (sphere-reject)
    
    % returns:
        % - feat: feature locations
        % - desc: feature descriptors
        
        
    % define histogram layout (bins) in spherical coordinates
    % it should be NUM_PHI = 2*NUM_THETA
    NUM_R = 2; % 2
    NUM_THETA = 3; % 3
    NUM_PHI = 6; % 6


       
    % preallocate space for features and descriptors    
    desc = nan(size(sample_pts, 1), NUM_R * NUM_PHI * NUM_THETA);
    feat = nan(size(sample_pts, 1), 3);
    
    tic
    parfor i = 1:size(sample_pts, 1) % PARFOR
        c = sample_pts(i, :);
        
        % return local points
        [pts_local, dists] = getLocalPoints(pts, R, c, min_pts);

        if ~ isempty(pts_local) 
            num_points = size(pts_local, 1);
            
            % ---- rejection based on PCA variances (== eigenvalues of
            % covariance matrix!!) of k nearest neighbors  
            if strcmp(K, 'all')
                k = num_points;
            else
                k = K;
                % sort points by distance to center for KNN
                [~, I] = sort(dists);
                pts_local = pts_local(I, :);
            end
            pts_k = pts_local(1:k, :);
            if ~(sum(thVar == 1) == 2) || ALIGN_POINTS
                [coeff, pts_lrf, variances] = pca(pts_k, 'Algorithm', 'eig'); 
                % continue if constraints are not met
                if (variances(1) / variances(2) < thVar(1)) || ...
                        (variances(2) / variances(3) < thVar(2))
                    continue
                end
            end

            % ---- use sign disambiguition method for aligned points
            
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
                pts_local = pts_local*coeff_unambig;
                %r = rand(1, 3)*2*pi;
                %pts_local = pts_local * eul2rotm(r);
                %pts_local = pts_local*coeff;
            end
            
            % ---- calculate the descriptor
            % center to mean/median point coordinates
            M1_L1 = mean(pts_local, 1)
            if CENTER
                pts_local = pts_local - M1_L1;
                c = c - M1_L1;
            end
            
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
            
            % flatten 3D histogram to 1D and normalize
            new_entry = reshape(counts, [], 1) / size(pts_local, 1);
            
            desc(i, :) = new_entry;
            feat(i, :) = c;
        end
    end
    
    % remove nan rows from desc and feat and return
    mask = find(~isnan(desc(:, 1))); % row indices of filled rows
    desc = desc(mask, :);
    feat = feat(mask, :);
    toc
end
