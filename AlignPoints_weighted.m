function [pts_aligned, coeff_unambig, c] = AlignPoints_weighted(pts)
%% AlignPoints_KNN.m
% PCA for alignment on the k (in %) nearest neighbors of the centroid
% returns the pca-aligned points (pts_aligned),
% the rotation matrix that was used (coeff_unambig)
% and the centroid 

    % --- 1) find median (L1) / mean (L2)
    c = mean(pts, 1);
    
    % --- 2) calculate weighted covariance matrix
    pts_rel = pts-c;
    dist_from_c = vecnorm(pts_rel, 2, 2);
    
    % get weight vector for points   
    R = 3.5;
    reverse_dist_from_c = R - dist_from_c;
    reverse_dist_from_c(reverse_dist_from_c < 0) = 0;
    %weighted_dist_from_c =  reverse_dist_from_c/ sum(reverse_dist_from_c);
    
    M = (reverse_dist_from_c .* pts_rel)' * pts_rel;
    
    % --- 3) get eigenvectors
    [coeff,~] = eig(M);   
    
    % ---- use sign disambiguition method for aligned points
    % centered points in local reference frame
    pts_lrf = pts_rel*coeff;
    
    k = size(pts, 1);
    
    % count number of points with positive sign and see if they
    % dominate 
    x_sign = sum(sign(pts_lrf(:, 1)) == 1) >= k/2;
    z_sign = sum(sign(pts_lrf(:, 3)) == 1) >= k/2;

    % map from {0, 1} to {-1, 1}
    x_sign = x_sign*2-1;
    z_sign = z_sign*2-1;

    %  get y sign so that rotation is proper
    y_sign = det(coeff .* [x_sign, 1, z_sign]);

    % apply signs to transform matrix (pca coefficients)
    coeff_unambig = coeff .* [x_sign, y_sign, z_sign];

    % --- 4) transform all points into the new coordinate system
    pts_aligned = pts*coeff_unambig;
end