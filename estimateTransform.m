%% function that estimates the 3D rigid transformation from at least 3 pts
function T = estimateTransform(pts1, pts2)
% this function estimates the rigid transformation T between pts1 and pts2
% based on PCA / covariance-eigenvectors. At least three points are needed for the 3D transformation. 

% returns T, such that [pts1, 1] * T = [pts2, 1]
    
    N = size(pts1, 1);
    
    C1 = sum(pts1, 1)./N;
    C2 = sum(pts2, 1)./N;
    
    M1 = pca(pts1, 'Algorithm', 'svd');
    M2 = pca(pts2, 'Algorithm', 'svd');
    
    % create third vector as cross-product, since pca returns 3x2 matrix
    % for only 3 points
    if N == 3
        M1(:, 3) = cross(M1(:, 1), M1(:, 2));
        M2(:, 3) = cross(M2(:, 1), M2(:, 2));
    end
    
    % center and align points with pca axes
    pts1_lrf = (pts1 - C1) * M1;
    pts2_lrf = (pts2 - C2) * M2;
    
    % now do sign disambiguition for both M1 and M2
    x_sign1 = sum(sign(pts1_lrf(:, 1)) == 1) >= N/2;
    z_sign1 = sum(sign(pts1_lrf(:, 3)) == 1) >= N/2;
    
    x_sign2 = sum(sign(pts2_lrf(:, 1)) == 1) >= N/2;
    z_sign2 = sum(sign(pts2_lrf(:, 3)) == 1) >= N/2;

    % map from {0, 1} to {-1, 1}
    x_sign1 = x_sign1*2-1;
    z_sign1 = z_sign1*2-1;
    
    x_sign2 = x_sign2*2-1;
    z_sign2 = z_sign2*2-1;

    %  get y sign so that M1 and M2 are actually rotation matrices
    y_sign1 = det(M1 .* [x_sign1, 1, z_sign1]);
    
    y_sign2 = det(M2 .* [x_sign2, 1, z_sign2]);
    
    %y_sign1 = M1(2, :) / cross(M1(1, :), M1(3, :)) * x_sign1 * z_sign1;
    
    %y_sign2 = M2(2, :) / cross(M2(1, :), M2(3, :)) * x_sign2 * z_sign2;

    % apply signs to transform matrix (pca coefficients)
    M1_unambig = M1 .* [x_sign1, y_sign1, z_sign1];
    
    M2_unambig = M2 .* [x_sign2, y_sign2, z_sign2];
    
    R = M1_unambig*M2_unambig';
    t = C2 - C1*R;
    
    T = eye(4);
    T(1:3, 1:3) = R;
    T(4, 1:3) = t;
end