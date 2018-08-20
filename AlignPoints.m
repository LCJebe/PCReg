function [pts_aligned, coeff_unambig] = AlignPoints(pts)

    [coeff, pts_lrf, ~] = pca(pts, 'Algorithm', 'eig'); 

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

    % transform all points into the new coordinate system
    pts_aligned = pts*coeff_unambig;
end