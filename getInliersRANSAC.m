%% run GetSphericalDescriptors first to load variables into workspace

DEBUG = false;

% matching points from model (loc1M) and surface(loc1S)
% to test with bad crop, use loc2M and loc2S

pts1 = loc2M; % Model
pts2 = loc2S; % Surface

%% debug: transform points to alignment manually and run RANSAC 
% (--> should return identity transform ideally)
if DEBUG
    load('Tform_align2.mat');
    pts2_postalign = [pts2, ones(size(pts2, 1), 1)]*Tform_align2;
    pts2_postalign = pts2_postalign(:, 1:3);
end

%% set RANSAC Parameters

coeff.minPtNum = 3;
coeff.iterNum = 3e4;
coeff.thDist = 1.0;
coeff.thInlrRatio = 0.1;

%% Perform RANSAC with rigid transform T and distance function
REFINE = true;
tic
if ~ DEBUG
    [T, inlierPtIdx] = ransac(pts1,pts2,coeff,@estimateTransform,@calcDists, REFINE);
else
    [T, inlierPtIdx] = ransac(pts1,pts2_postalign,coeff,@estimateTransform,@calcDists, REFINE);
end
toc


%% Use returend T to align pts1 (Model) with pts2 (Surface)
pts1_aligned = [pts1, ones(size(pts1, 1), 1)] * T;
pts1_aligned = pts1_aligned(:, 1:3);


%% DEBUG: If manually aligned post-matching, use simple distance check for
% inliers
if DEBUG
    thDist = 1;
    d1 = vecnorm(pts2_postalign - pts1, 2, 2);
    inliers_postalign = length(find(d1 < thDist));
end

%% function that calculates the distance between points after transform T
function d = calcDists(T,pts1,pts2)
    %	Project PTS1 to PTS1_trans using the rigid transform T, then calcultate the distances between
    %	PTS2 and PTS1_trans

    pts1(:, 4) = 1;
    pts1_trans = pts1*T;
    pts1_trans = pts1_trans(:, 1:3);
    d = sum((pts2-pts1_trans).^2,2);
end