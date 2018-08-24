function [pts_aligned, coeff_unambig, c] = AlignPoints_KNN(pts)
    % returns the pca-aligned points (pts_aligned)
    % and the rotation matrix that was used (coeff_unambig)

    % --- 1) find median (L1) / mean (L2)
    c = mean(pts, 1);
    
    % --- 2) get the k nearest neighbors from the centroid
    k = 0.85; % fraction of points considered. (rest assumed outlier)
    K = round(size(pts,1) * k);
    pts_rel = pts - c;
    dists = vecnorm(pts_rel, 2, 2);
    [~, I] = sort(dists);
    pts_sorted = pts_rel(I, :);
    pts_k = pts_sorted(1:K, :);
    
    
    % --- 3) use those points for alignment (get transform)
    [coeff, pts_lrf, ~] = pca(pts_k, 'Algorithm', 'eig'); 

    % ---- use sign disambiguition method for aligned points
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