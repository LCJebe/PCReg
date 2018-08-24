function [pts_aligned, coeff_unambig, c] = AlignPoints_v2(pts)
    % returns the pca-aligned points (pts_aligned)
    % and the rotation matrix that was used (coeff_unambig)

    % --- 1) find median (L1) / mean (L2)
    c = mean(pts, 1);
    
    % --- 2) get local neighborhood around this point with radius r and
    % make sure there is at least a bare minimum of points included
    r = 2.0;
    [pts_rel, ~] = getLocalPoints(pts, r, c, 25, inf);
    
    if isempty(pts_rel)
        pts_aligned = [];
        coeff_unambig = [];
    else
        % --- 3) use those points for alignment (get transform)
        [coeff, pts_lrf, ~] = pca(pts_rel, 'Algorithm', 'eig'); 

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
end