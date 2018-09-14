function pts_tf = quickTF(pts, TF)
% transforms a set of 3D points with a given 4x4 transform matrix
% points are of shape N x 3

    pts = [pts, ones(size(pts, 1), 1)];
    pts_tf = pts * TF;
    pts_tf = pts_tf(:, 1:3);
end

