close all

%% load model and aligned surface
path = 'Data/PointClouds/';

pcSurface = pcread(strcat(path, 'Surface_DS2_alignedM.pcd'));
pcModel = pcread(strcat(path, 'GoodCropSmoothUp3_large.pcd'));

%% options
ALIGN = true;
MAX_MATCHES = 8;

%% sample points (same for surface and model)
d = 1;
sample_pts = pcRandomUniformSamples(pcSurface, d);

%% select points, where surface and model are not empty
R = 3.5;
min_pts = 501;
max_pts = 8*min_pts;

valid_pts = [];
count_validS = 0;
count_validM = 0;
for i = 1:size(sample_pts ,1)
    c = sample_pts(i, :);
    % for Surface
    [local_ptsS, ~] = getLocalPoints(pcSurface.Location, R, c, min_pts, max_pts);
    % for Model
    [local_ptsM, ~] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);
    
    % only keep when both are valid (not empty)
    if (~isempty(local_ptsS)) && (~isempty(local_ptsM))
        valid_pts = [valid_pts; c];
    elseif ~isempty(local_ptsS)
        count_validS = count_validS + 1;
    elseif ~isempty(local_ptsM)
        count_validM = count_validM + 1;
    end
end
clear local_ptsS local_ptsM

% print out some intermediate statistics
num_samples = size(sample_pts ,1);
num_valid = size(valid_pts, 1);
perc_valid = double(num_valid) / num_samples * 100;
fprintf('Sampled %d points, %d are valid (%0.2f %%)\n', ...
   num_samples, num_valid, perc_valid);


%% now the fun part: visualize ground truth matches

% for better display: generate the 8 points that span the cube defined by R
add_pts = R*[-1,-1,-1;...
             -1,-1, 1;...
             -1, 1,-1;...
             -1, 1, 1;...
              1,-1,-1;...
              1,-1, 1;...
              1, 1,-1;...
              1, 1, 1];

% figure settings
screensize = get( 0, 'Screensize' );
if ~ALIGN
    figpos = [screensize(3)/3, 75, screensize(3)/3, screensize(4)-150];
else 
    figpos = [screensize(3)/6, 75, 2*screensize(3)/3, screensize(4)-150];
end

num_display = min([MAX_MATCHES, num_valid]);

for i = 1:num_display
    
    % get a random selection of the matches
    matches_vis = randperm(size(valid_pts, 1),num_display);
    c = valid_pts(matches_vis(i), :);
    
    % get points
    [ptsMatchSurface, ~] = getLocalPoints(pcSurface.Location, R, c, min_pts, max_pts);
    [ptsMatchModel, ~] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);
    
    % optional: align to local reference frame
    if ALIGN
        [ptsMatchSurface_aligned, ~] = AlignPoints(ptsMatchSurface);
        [ptsMatchModel_aligned, ~] = AlignPoints(ptsMatchModel);
    end
    
    % add frame points for display
    ptsMatchSurface = [ptsMatchSurface; add_pts];
    ptsMatchModel = [ptsMatchModel; add_pts];
    if ALIGN
        ptsMatchSurface_aligned = [ptsMatchSurface_aligned; add_pts];
        ptsMatchModel_aligned = [ptsMatchModel_aligned; add_pts];
    end
    
    % visualize
    if ~ ALIGN
        row = mod(i-1, 4)+1;
        if row ==1
           fig_h = figure();
           set(fig_h,'Position',figpos)
        end
        subplot(4, 2, 2*row-1);
        pcshow(ptsMatchSurface, 'MarkerSize', 50);
        if row == 1; title('Surface'); end
        subplot(4, 2, 2*row);
        pcshow(ptsMatchModel, 'MarkerSize', 50);
        if row == 1; title('Model'); end
    else
        row = mod(i-1, 4)+1;
        if row ==1
           fig_h = figure();
           set(fig_h,'Position',figpos)
        end
        subplot(4, 4, 4*row-3);
        pcshow(ptsMatchSurface, 'MarkerSize', 50);
        if row == 1; title('Surface'); end
        subplot(4, 4, 4*row-2);
        pcshow(ptsMatchModel, 'MarkerSize', 50);
        if row == 1; title('Model'); end
        subplot(4, 4, 4*row-1);
        pcshow(ptsMatchSurface_aligned, 'MarkerSize', 50);
        if row == 1; title('Surface aligned'); end
        subplot(4, 4, 4*row);
        pcshow(ptsMatchModel_aligned, 'MarkerSize', 50);
        if row == 1; title('Model aligned'); end
    end
end




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