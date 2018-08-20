pcModel = pcread(strcat(path, 'GoodCropSmoothUp3_large.pcd'));

%% sample points
d = 1;
sample_pts = pcRandomUniformSamples(pcSurface, d);

%% select points, where surface and model are not empty
R = 3.5;
min_pts = 501;
max_pts = 8*min_pts;

% run get local points and time
tic
for i = 1:size(sample_pts ,1)
    c = sample_pts(i, :);
    [local_ptsM, ~] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);
end
toc 

tic
for i = 1:size(sample_pts ,1)
    c = sample_pts(i, :);
    [local_ptsM, ~] = getLocalPoints_v2(pcModel.Location, R, c, min_pts, max_pts);
end
toc 


%% helper function: random-uniformly sample point cloud
function sample_pts = pcRandomUniformSamples(pcIn, d)
    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1);
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1);
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1);
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers
    sample_pts = rand(num_pts, 3);
    
    % scale numbers so that they fit into the correct range
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)];
end