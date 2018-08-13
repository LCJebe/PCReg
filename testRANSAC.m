% script to test implementation of RANSAC
clear all

%% load test pointcloud (teapot)

ptCloud = pcread('teapot.ply');
pts = ptCloud.Location;
num_pts = size(pts, 1);

%% define transformation and transform points

% specify Rotation and Translation
r = [1.5,-1.2,0.8];
t = [1,2,3];

% rotate and translate points accordingly
R = eul2rotm(r);
pts_tf = pts*R + t;


%% get correct transformation matrix

% estimate transformation between pts and pts_tf
T_true1 = estimateTransform(pts, pts_tf);

% hardcode correct matrix
T_true2 = eye(4);
T_true2(1:3, 1:3) = R;
T_true2(4, 1:3) = t;

% use this transformation to transform pts to pts_tf
pts_tf1 = [pts, ones(num_pts, 1)]*T_true1;
pts_tf2 = [pts, ones(num_pts, 1)]*T_true2;

% check the error
error1 = mean(sqrt(sum((pts_tf-pts_tf1(:, 1:3)).^2, 2)));
error2 = mean(sqrt(sum((pts_tf-pts_tf2(:, 1:3)).^2, 2)));

% hardcode correct backtransformation matrix as well
T_back = eye(4);
T_back(1:3, 1:3) = R';
T_back(4, 1:3) = -t*R';

% transform pts_tf2 back to pts
pts_restored = pts_tf2*T_back;

% check the error
error3 = mean(sqrt(sum((pts-pts_restored(:, 1:3)).^2, 2)));


%% optional: introduce gaussian noise to points
sigma = 0.1;
pts_noisy = pts + randn(size(pts))*sigma;


%% use RANSAC to find correct transformation
loc1M = pts_noisy;
loc1S = pts_tf2(:, 1:3);

% now call getInliersRANSAC.m !!