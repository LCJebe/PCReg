%% run GetSphericalDescriptors first to load everything into the workspace

%% define parameter sweep
start = 0.3;
until = 0.7;
steps = 30;
maxRatio = start:(until-start)/steps:until;
thMatch = logspace(-1, 0, 31);


%% initialize result metrics

precision = zeros(length(maxRatio), length(thMatch));
goodInliers = zeros(length(maxRatio), length(thMatch));
badInliers = zeros(length(maxRatio), length(thMatch));

%% %% Match features between Surface and Model / Random Crop
tic
for i = 1:length(maxRatio)
    for j = 1:length(thMatch)

        matchesModel = matchFeatures(descSurfaceW, descModelW, ...
                'Method', 'Approximate', ...
                'MatchThreshold', thMatch(j), ... 
                'MaxRatio', maxRatio(i), ... 
                'Metric', 'SSD', ...
                'Unique', true); 

        matchesRand = matchFeatures(descSurfaceW, descRandW, ...
                'Method', 'Approximate', ...
                'MatchThreshold', thMatch(j), ... 
                'MaxRatio', maxRatio(i), ... 
                'Metric', 'SSD', ...
                'Unique', true); 

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
        precision(i, j) = inliersPrecision;
        goodInliers(i, j) = inliers1;
        badInliers(i, j) = inliers2;
    end
    toc
end

%% Plots
close all

[X, Y] = meshgrid(thMatch, maxRatio);

figure()
surf(X, Y, precision, 'FaceColor', 'interp');
title("Precision (%)");
xlabel("Nearest Neighbor Ratio");
ylabel("Distance Threshold");
colorbar
set(gca, 'XScale','log');


figure()
surf(X, Y, goodInliers, 'FaceColor', 'interp');
title("Good Inliers (#)");
xlabel("Nearest Neighbor Ratio");
ylabel("Distance Threshold");
colorbar
set(gca, 'XScale','log');

figure()
surf(X, Y, goodInliers.*precision/100, 'FaceColor', 'interp');
title("Good Inliers (#) x Precision");
xlabel("Nearest Neighbor Ratio");
ylabel("Distance Threshold");
colorbar
set(gca, 'XScale','log');