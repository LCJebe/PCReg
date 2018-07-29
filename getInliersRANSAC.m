%% run GetSphericalDescriptors first to load variables into workspace

% matching points from model (loc1M) and surface(loc1S);
pts1 = loc1M;
pts2 = loc1S;

%% set RANSAC Parameters

coeff.minPtNum = 3;
coeff.iterNum = 2e4;
coeff.thDist = 2;
coeff.thInlrRatio = .05;

%% Perform RANSAC with rigid transform T and distance function
num_trials = 10;

tic
[T, inlierPtIdxm] = ransac(pts1,pts2,coeff,@estimateTransform,@calcDists);
toc

%% function that calculates the distance between points after transform T
function d = calcDists(T,pts1,pts2)
    %	Project PTS1 to PTS1_trans using the rigid transform T, then calcultate the distances between
    %	PTS2 and PTS1_trans

    pts1(:, 4) = 1;
    pts1_trans = pts1*T;
    pts1_trans = pts1_trans(:, 1:3);
    d = sum((pts2-pts1_trans).^2,2);
    
end

