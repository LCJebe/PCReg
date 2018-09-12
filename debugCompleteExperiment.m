%% load data    
load('Data/Descriptors/featSurface0.3.mat');
load('Data/Descriptors/descSurface0.3.mat'); 
load('Data/Descriptors/featModel0.3.mat');
load('Data/Descriptors/descModel0.3.mat'); 

centerSurface = [34.8970   21.9820   25.4546];
R_desc = 7.0882;

%% match the ideal sphere
% matching options
par.UNNORMALIZE = true;
par.norm_factor = 2; % 2

par.CHANGE_METRIC = true;
par.metric_factor = 0.6; % 0.45 / 0.6

par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold = 10; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.99; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SAD'; % SSD (default) for L2, SAD for L1
par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

par.VERBOSE = 1;    

% only one sphere around the center
c = centerSurface + [0 0 0];

good_samples = ...
   [35.5690   21.7328   21.5702;
   31.5690   21.7328   23.5702;
   33.5690   19.7328   23.5702;
   33.5690   21.7328   23.5702;
   33.5690   23.7328   23.5702;
   35.5690   19.7328   23.5702;
   35.5690   21.7328   23.5702;
   35.5690   23.7328   23.5702;
   37.5690   19.7328   23.5702;
   37.5690   21.7328   23.5702;
   37.5690   23.7328   23.5702;
   31.5690   21.7328   25.5702;
   31.5690   23.7328   25.5702;
   33.5690   19.7328   25.5702;
   33.5690   21.7328   25.5702;
   33.5690   23.7328   25.5702;
   33.5690   25.7328   25.5702;
   35.5690   19.7328   25.5702;
   35.5690   21.7328   25.5702;
   35.5690   23.7328   25.5702;
   35.5690   25.7328   25.5702;
   37.5690   19.7328   25.5702;
   37.5690   21.7328   25.5702;
   37.5690   23.7328   25.5702;
   31.5690   21.7328   27.5702;
   33.5690   19.7328   27.5702;
   33.5690   21.7328   27.5702;
   33.5690   23.7328   27.5702;
   35.5690   19.7328   27.5702;
   35.5690   21.7328   27.5702;
   35.5690   23.7328   27.5702;
   37.5690   21.7328   27.5702;
   37.5690   23.7328   27.5702];

%% RANSAC options
options.minPtNum = 3; 
options.iterNum = 1e4; 
options.thDist = 0.2; 
options.thInlrRatio = 0.1; 
options.REFINE = true;
options.VERBOSE = 1;

% mask of descriptors that lie within sphere
good_mask = false(size(good_samples, 1), 1);
num_matches = zeros(size(good_samples, 1), 1);
for i = 1:size(good_samples, 1)
    c = good_samples(i, :);
    
    mask = getDescriptorIndices(featModel, c, R_desc, 0);

    % based on the mask, return the relevant sample points and descriptors
    maskF = cat(2, mask, mask, mask);
    featCur = reshape(featModel(maskF), [], 3);  
    maskD = repmat(mask, 1, size(descModel, 2));
    descCur = reshape(descModel(maskD), [], size(descModel, 2));

    % now match with the specified region
    matchesModel = getMatches(descSurface, descCur, par);

    fprintf('Found %d putative matches\n', size(matchesModel, 1));


    % run ransac
    pts1 = featSurface(matchesModel(:, 1), :);
    pts2 = featCur(matchesModel(:, 2), :);

    [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = ...
                ransac(pts1,pts2,options,@estimateTransform,@calcDists);
            
    % record sucesses and display sphere centers (Surface center in red,
    % successes in green, unsuccessful trials in black)
    if numSuccess > 0
        good_mask(i) = 1;
    end
    
    % record some statistics
    num_matches(i) = size(pts1, 1);

end

%% Display results for debugging (which spheres succeed?)
% record sucesses and display sphere centers (Surface center in red,
% successes in green, unsuccessful trials in black)

goodLocations = reshape(good_samples(repmat(logical(good_mask), 1, 3)), [], 3);
goodColor = [0, 255, 0] .* ones(size(goodLocations));
badLocations = reshape(good_samples(repmat(~good_mask, 1, 3)), [], 3);
badColor = [0, 0, 0] .* ones(size(badLocations));

color = [goodColor; badColor; [255, 0, 0]];
points = [goodLocations; badLocations; centerSurface];

pcshow(pointCloud(points, 'Color', color), 'MarkerSize', 50);

%% helper function: descriptors within certain distance of given center point
function mask = getDescriptorIndices(featModel, current_center, R_desc, margin)
    feat_rel = featModel - current_center;
    dists = vecnorm(feat_rel, 2, 2);
    mask = dists < (R_desc + margin); % note that margin < 0
end
%% function that calculates the distance between points after transform T
% used as distance metric for RANSAC
function d = calcDists(T,pts1,pts2)
    %	Project PTS2 to PTS2_trans using the rigid transform T, then calcultate the distances between
    %	PTS1 and PTS2_trans
    pts2(:, 4) = 1;
    pts2_trans = pts2*T;
    pts2_trans = pts2_trans(:, 1:3);
    d = sum((pts1-pts2_trans).^2,2);
end