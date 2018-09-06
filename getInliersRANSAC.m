%% run GetSphericalDescriptors first to load variables into workspace!

%% getInliersRANSAC
% After matching of descriptors has been performed, RANSAC (using a 3D
% rigid body transformation) is an easy way to obtain the inlier
% (==correct) matches. 


% select loc1M / loc1S to match the query with the correct crop
% select loc2M / loc2S to match the query with the wrong crop. Ideally,
% RANSAC should fail in this case and not return a transformation
pts1 = loc1M; % Model
pts2 = loc1S; % Surface

%% set RANSAC Parameters
% number of sampled points. minimum is 3 for rigid body transformation
coeff.minPtNum = 3; 

% number of iterations
coeff.iterNum = 2e4; 

% distance (in world units) below which matches are considered inliers
coeff.thDist = 0.5; 

% percentage of matches that are inliers needed to call the transformation a success
coeff.thInlrRatio = 0.1; 


%% Perform RANSAC with rigid transform T and distance function
% REFINE: find the transformation again using all inliers, if successful
coeff.REFINE = true;

tic
[T, inlierPtIdx] = ransac(pts1,pts2,coeff,@estimateTransform,@calcDists);
toc

T
%% Use returend T to align pts1 (Model) with pts2 (Surface)
if ~isempty(T)
    pts1_aligned = [pts1, ones(size(pts1, 1), 1)] * T;
    pts1_aligned = pts1_aligned(:, 1:3);
end


%% function that calculates the distance between points after transform T
function d = calcDists(T,pts1,pts2)
    %	Project PTS2 to PTS2_trans using the rigid transform T, then calcultate the distances between
    %	PTS1 and PTS2_trans

    pts2(:, 4) = 1;
    pts2_trans = pts2*T;
    pts2_trans = pts2_trans(:, 1:3);
    d = sum((pts1-pts2_trans).^2,2);
end