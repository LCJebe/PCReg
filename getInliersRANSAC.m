%% run GetSphericalDescriptors first to load variables into workspace

% matching points from model (loc1M) and surface(loc1S)
% to test with bad crop, use loc2M and loc2S

pts1 = loc1M; % Model
pts2 = loc1S; % Surface

%% set RANSAC Parameters

coeff.minPtNum = 3;
coeff.iterNum = 4e4;
coeff.thDist = 0.5;
coeff.thInlrRatio = 0.1;

%% Perform RANSAC with rigid transform T and distance function
REFINE = true;
tic
[T, inlierPtIdx] = ransac(pts1,pts2,coeff,@estimateTransform,@calcDists, REFINE);
toc

T
%% Use returend T to align pts1 (Model) with pts2 (Surface)
pts1_aligned = [pts1, ones(size(pts1, 1), 1)] * T;
pts1_aligned = pts1_aligned(:, 1:3);

%% function that calculates the distance between points after transform T
function d = calcDists(T,pts1,pts2)
    %	Project PTS1 to PTS1_trans using the rigid transform T, then calcultate the distances between
    %	PTS2 and PTS1_trans

    pts1(:, 4) = 1;
    pts1_trans = pts1*T;
    pts1_trans = pts1_trans(:, 1:3);
    d = sum((pts2-pts1_trans).^2,2);
end