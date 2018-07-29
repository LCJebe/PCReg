%% run GetSphericalDescriptors first to load everything into the workspace

%% define parameter sweep
thMatch = logspace(-1, 1, 50);


%% initialize result metrics

precision = zeros(length(thMatch), 1);
goodInliers = zeros(length(thMatch), 1);
badInliers = zeros(length(thMatch), 1);

%% %% Match features between Surface and Model / Random Crop

for i = 1:length(thMatch)

    % Define matching algorithm parameters
    par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
    par.MatchThreshold =  thMatch(i); % 1.0 (default) Percent Value (0 - 100) for distance-reject
    par.MaxRatio = 0.8; % 0.6 (default) nearest neighbor ambiguity rejection
    par.Metric =  'SSD'; % SSD (default) for L2, SAD for L1
    par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

    matchesModel = matchFeatures(descSurfaceW, descModelW, ...
            'Method', par.Method, ...
            'MatchThreshold', par.MatchThreshold, ... 
            'MaxRatio', par.MaxRatio, ... 
            'Metric', par.Metric, ...
            'Unique', par.Unique); 

    matchesRand = matchFeatures(descSurfaceW, descRandW, ...
            'Method', par.Method, ...
            'MatchThreshold', par.MatchThreshold, ... 
            'MaxRatio', par.MaxRatio, ... 
            'Metric', par.Metric, ...
            'Unique', par.Unique); 

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
    loc2M = featRand(matchesRand(:, 2), :);
    loc2S = featSurface(matchesRand(:, 1), :);
    d2 = vecnorm(loc2M - loc2S, 2, 2);
    inliers2 = length(find(d2 <= maxDist));
    
    % save results for this run
    precision(i) = inliersPrecision;
    goodInliers(i) = inliers1;
    badInliers(i) = inliers2;
end

%% Plots
close all

figure()
semilogx(thMatch, precision, '-*');
%title("Precision (%)");
xlabel("Match Threshold");
grid;

%figure()
hold on
semilogx(thMatch, goodInliers, '-*');
%title("Good Inliers (#)");
xlabel("Match Threshold");
%grid;
legend("Precision (%)", "Good Inliers (#)");

figure()
semilogx(thMatch, badInliers);
title("Bad Inliers");
xlabel("Match Threshold");
grid;