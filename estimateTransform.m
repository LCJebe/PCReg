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
    
    R = M1*M2';
    t = C2 - C1*R;
    
    T = eye(4);
    T(1:3, 1:3) = R;
    T(4, 1:3) = t;
end