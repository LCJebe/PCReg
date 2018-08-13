%% Load variables into workspace first by calling GetSphericalDescriptors!s

%% Sample combintions of m2, m3 on 1/4 circle

steps = 100;
theta = 0:pi/2/steps:pi/2;
theta_deg = rad2deg(theta);


%%
weights.M1 = 0;
weights.M3 = 0;

precision = zeros(length(theta), 1);
goodInliers = zeros(length(theta), 1);

tic
for t = 1:length(theta)

    weights.M2 = cos(theta(t));
    weights.M4 = sin(theta(t)); 

    descSurfaceW = applyMomentWeight(descSurface, weights);
    descModelW = applyMomentWeight(descModel, weights);
    descRandW = applyMomentWeight(descRand, weights);

    %% Match features between Surface and Model / Random Crop

    % Define matching algorithm parameters
    par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
    par.MatchThreshold =  1; % 1.0 (default) Percent Value (0 - 100) for distance-reject
    par.MaxRatio = 0.8; % 0.6 (default) nearest neighbor ambiguity rejection
    par.Metric =  'SSD'; % SSD (default) for L2, SAD for L1
    par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

    matchesModel = matchFeatures(descSurfaceW, descModelW, ...
            'Method', par.Method, ...
            'MatchThreshold', par.MatchThreshold, ... 
            'MaxRatio', par.MaxRatio, ... 
            'Metric', par.Metric, ...
            'Unique', par.Unique); 

%     matchesRand = matchFeatures(descSurfaceW, descRandW, ...
%             'Method', par.Method, ...
%             'MatchThreshold', par.MatchThreshold, ... 
%             'MaxRatio', par.MaxRatio, ... 
%             'Metric', par.Metric, ...
%             'Unique', par.Unique); 

    %% Get distance between matching points
    % this makes sense, because the pointclouds are already aligned. The
    % distance between matching INLIERS of the model will thus be small and a
    % simple distance threshold can be used to determine whether a match is an
    % inlier 

    % maxDist specifies matching distance to count inliers 
    maxDist = 1; 

    % matches of surface and model
    loc1M = featModel(matchesModel(:, 2), :);
    loc1S = featSurface(matchesModel(:, 1), :);
    d1 = vecnorm(loc1M - loc1S, 2, 2);
    inliers1 = length(find(d1 < maxDist));

    % precision
    inliersPrecision = inliers1/size(d1, 1)*100;

    % matches of surface and random crop
%     loc2M = featRand(matchesRand(:, 2), :);
%     loc2S = featSurface(matchesRand(:, 1), :);
%     d2 = vecnorm(loc2M - loc2S, 2, 2);
%     inliers2 = length(find(d2 <= maxDist));

    % fill results for curve
    precision(t) = inliersPrecision;
    goodInliers(t) = inliers1;

end 

toc
%% Plots
close all;

figure()
plot(theta_deg, precision)
title("Precision m2 --- m4 tradeoff")
xlabel("Theta")
ylabel("Precision")
grid;

figure()
plot(theta_deg, goodInliers)
title("Good Inliers m2 --- m4 tradeoff")
xlabel("Theta")
ylabel("Num Inliers")
grid;

%% helper function that applies weight to descriptor
function descW = applyMomentWeight(desc, weights)
    descW = desc;
    descW(:, 1:6) = desc(:, 1:6)*weights.M1;
    descW(:, 7:12) = desc(:, 7:12)*weights.M2;
    descW(:, 13:22) = desc(:, 13:22)*weights.M3;
    descW(:, 23:31) = desc(:, 23:31)*weights.M4;
end