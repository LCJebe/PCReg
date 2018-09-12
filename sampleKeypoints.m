function sample_pts = sampleKeypoints(pcIn, parameters)
%% sample keypoints given a pointcloud and sampling parameters

% 'UNIFORM' for a regular grid (uniform sampling) for each pointcloud
% 'UNIFORM_SAME' for the same regular grid on all point clouds
% 'RANDOM_POINTS' for random points in point cloud
% 'RANDOM_UNIFORM' for random uniformly distributed points
% 'RANDOM_UNIFORM_SPHERICAL' for random uniform points in sphere defined by center and radius
% 'ALL' for all points in point cloud

    % unpack parameters
    sampling_method = parameters.sampling_method;
    d = parameters.d;
    margin = parameters.margin;
    if isfield(parameters, 'sample_frac')
        sample_frac = parameters.sample_frac;
    end

    % chose function to execute
    if strcmp(sampling_method, 'UNIFORM')    
        sample_pts = pcUniformSamples(pcIn, d);

    elseif strcmp(sampling_method, 'RANDOM_UNIFORM')
        sample_pts = pcRandomUniformSamples(pcIn, d, margin);

    elseif strcmp(sampling_method, 'RANDOM_UNIFORM_SPHERICAL')    
        sample_pts = pcRandomUniformSphericalSamples(pcIn, d, margin);

    elseif strcmp(sampling_method, 'RANDOM_POINTS')  
        try
            sample_pts = pcRandomPoints(pcIn, sample_frac);
        catch
            error('For method RANDOM_POINTS please sepcify a fraction in the field parameters.sample_frac');
        end

    elseif strcmp(sampling_method, 'ALL')    
        sample_pts = pcSurface.Location;

    else
        error('Specified Method "%s" not yet implemented or does not exist', sampling_method);
    end

end



%%% --- define local functions --- %%%

%% helper function: uniformly sample point cloud
function sample_pts = pcUniformSamples(pcIn, d) 
    % create meshgrid of points	  
    x = pcIn.XLimits(1):d:pcIn.XLimits(2);
    y = pcIn.YLimits(1):d:pcIn.YLimits(2);
    z = pcIn.ZLimits(1):d:pcIn.ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
end

%% helper function: random-uniformly sample point cloud
function sample_pts = pcRandomUniformSamples(pcIn, d, margin)
    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1) + 2*margin;
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1) + 2*margin;
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1) + 2*margin;
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers in range [0, 1]
    sample_pts = rand(num_pts, 3); 
    
    % scale numbers so that they fit into the correct range
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)] - margin;
end

%% helper function: random-uniformly sample point cloud
function sample_pts = pcRandomUniformSphericalSamples(pcIn, d, margin)

    if margin > 0
        error('Margin is expected to be < 0 for this implementation of spherical sampling');
    end

    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1);
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1);
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1);
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers in range [0, 1]
    sample_pts = rand(num_pts, 3); 
    
    % scale numbers so that they fit into the correct range, without using
    % margin yet
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)];
    
    % now get center and only use points that are within "reach"
    center = [(pcIn.XLimits(2) + pcIn.XLimits(1))/2, ...
              (pcIn.YLimits(2) + pcIn.YLimits(1))/2, ...
              (pcIn.ZLimits(2) + pcIn.ZLimits(1))/2];
          
    pts_rel = sample_pts - center;
    dists = vecnorm(pts_rel, 2, 2);
    mask = dists < (rangeX/2 + margin); % note that margin < 0
    
    % based on the mask, return the relevant sample points
    mask = cat(2, mask, mask, mask);
    sample_pts = reshape(sample_pts(mask), [], 3);  
end

%% helper function: select random points from point cloud
function sample_pts = pcRandomPoints(pcIn, sample_frac)
    nPoints = size(pcIn.Location, 1);
    sIdx = randsample(nPoints,floor(nPoints*sample_frac));
    sample_pts = pcIn.Location(sIdx, :);
end