%% Load variables into workspace first by calling GetSphericalDescriptors!

%% Sample combintions of m2, m3, m4 on 1/8 sphere surface

% sample uniform theta
tsteps = 20;
theta = 0:pi/2/tsteps:pi/2;
theta_deg = rad2deg(theta);

% sample phi in arcos domain
psteps = 20;
punif = 0:1/psteps:1;
phi = acos(1-punif);
phi_deg = rad2deg(phi);

% grid for surf plot
[theta_grid, phi_grid] = meshgrid(theta_deg, phi_deg);

%%

precision = zeros(length(theta), length(phi));
goodInliers = zeros(length(theta), length(phi));

tic
for p = 1:length(phi)
    parfor t = 1:length(theta)

        weights = [0, sin(phi(p)) * cos(theta(t)), ...
            sin(phi(p)) * sin(theta(t)), ...
            cos(phi(p))];  

        descSurfaceW = applyMomentWeight(descSurface, weights);
        descModelW = applyMomentWeight(descModel, weights);

        %% Match features between Surface and Model / Random Crop

        matchesModel = matchFeatures(descSurfaceW, descModelW, ...
                'Method', 'Approximate', ...
                'MatchThreshold', 0.6, ... 
                'MaxRatio', 0.6, ... 
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

        % fill results for curve
        precision(t, p) = inliersPrecision;
        goodInliers(t, p) = inliers1;

    end 
    toc
end
%% Plots
close all;

figure()
surf(theta_grid, phi_grid, precision, 'FaceColor', 'interp')
title("Precision (%)")
xlabel("Theta")
ylabel("Phi")
colorbar

figure()
surf(theta_grid, phi_grid, goodInliers, 'FaceColor', 'interp')
title("# of Good Inliers")
xlabel("Theta")
ylabel("Phi")
colorbar

%% helper function that applies weight to descriptor
function descW = applyMomentWeight(desc, weights)
    descW = desc;
    descW(:, 1:6) = desc(:, 1:6)*weights(1);
    descW(:, 7:12) = desc(:, 7:12)*weights(2);
    descW(:, 13:22) = desc(:, 13:22)*weights(3);
    descW(:, 23:31) = desc(:, 23:31)*weights(4);
end