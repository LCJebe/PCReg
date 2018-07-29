% create three or four 3D points
pts = [1,5,7;
       4,9,3;
       9,3,4;
       1,2,4];

% specify Rotation and Translation
r = [0.1,0.2,0.3];
t = [1,2,3];

% rotate and translate points accordingly
R = eul2rotm(r);

pts_tf = pts*R + t;

% estimate transformation between pts and pts_tf
T = estimateTransform(pts_tf, pts);

% use transformation to transform pts_tf back to pts
pts_restored = [pts_tf, ones(4, 1)]*T;

% check the error
error = sqrt(sum((pts-pts_restored(:, 1:3)).^2))


