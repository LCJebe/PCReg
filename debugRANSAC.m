%% create N random points in 3D with random colors
N = 3;
pts = randn(N, 3);
pts_h = [pts, ones(N, 1)];

colors = uint8(rand(N, 3)*255);

% create a random rigid transformation T and its inverse Tinv
r = rand(1, 3)*2*pi;
t = randn(1, 3);
R = eul2rotm(r, 'XYZ');
T = eye(4);
T(1:3, 1:3) = R;
T(4, 1:3) = t;
Tinv = inv(T);

% transform points using the 3D rigid transformation
pts_tf_h = pts_h*T;
pts_tf = pts_tf_h(:, 1:3);

% estimate transformation using estimateTransform.m and get estimation error
T_est_ideal = estimateTransform(pts, pts_tf);

pts_back_tf_ideal_h = pts_tf_h * T_est_ideal;
pts_back_tf_ideal = pts_back_tf_ideal_h(:, 1:3);

% add noise to transformed points with standard deviation std
std = 0.1;
noise = std*randn(N, 3);
pts_tf_noisy = pts_tf + noise;

% estimate transform under noise
T_est_noisy = estimateTransform(pts, pts_tf_noisy);

pts_back_tf_noisy_h = pts_tf_h * T_est_noisy;
pts_back_tf_noisy = pts_back_tf_noisy_h(:, 1:3);

est_error = norm(T*T_est_noisy - eye(4), 'fro');

DISP = false;
if DISP
    pcshow(pointCloud([pts; pts_back_tf_noisy], 'Color', [colors; colors]), 'MarkerSize', 1000);
end

%% define parameters for ransac
coeff.minPtNum = 3; 
coeff.iterNum = 1e3;
coeff.thDist = 0.1; 
coeff.thInlrRatio = 0.5; 

%% run ransac to find correct transformation
REFINE = true;

[T_est, inlierPtIdx] = ransac(pts,pts_tf,coeff,@estimateTransform,@calcDists, REFINE);


%% function that calculates the distance between points after transform T
function d = calcDists(T,pts1,pts2)
    %	Project PTS2 to PTS2_trans using the rigid transform T, then calcultate the distances between
    %	PTS1 and PTS2_trans

    pts2(:, 4) = 1;
    pts2_trans = pts2*T;
    pts2_trans = pts2_trans(:, 1:3);
    d = sum((pts1-pts2_trans).^2,2);
end