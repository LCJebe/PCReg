%% Read stereo image input
%clearvars;
%close all;
l = imread('Data/Stereo/Nah2/02_L.tif');
r = imread('Data/Stereo/Nah2/02_R.tif');

%% Define METHOD
% either 'SPARSE' or 'DENSE'
% 'SPARSE' Method useds keypoint detection and matching to obtain a
% depthmap
% 'DENSE' Method uses Block matching provided my MATLAB. Recommended. 
METHOD = 'DENSE';

%% Define Camera Parameters and Information

% Focal Length in [pixels], the same in all directions for both cameras
FocalLength = 47805.168; 

% Rotation of camera 2 with respect to camera 1 in [rad]
theta_x = 0;
theta_y = 0.087; % 4.98 deg
theta_z = 0;

% Translation of camera 1 with respect to camera 2
T = [35, 0, 0];
 
% Radial Distortion Coefficients for camera 1 and camera 2
% Follow the Formula: sum(i=1 -> inf) (k_i r^(2i+1))
% Approximation with polynomial:
% k1 * r^3 + k2 * r^5, where r is the NORMALIZED distance from the center
% OR: OpenCV:  k1 * r^2 + k2 * r^4, where r is the NORMALIZED distance from the center
k1_left = 0.00346416;
k2_left =-0.00078215;

k1_right = 0.00168054;
k2_right = 0.00042590;

% pixel width in world units [mm]:
PixelWidth = 8.25e-3;

%% Define StereoParameters Object from given Information
FL = [FocalLength, FocalLength];
IS = [1080, 1920];
PP = IS/2;

% radial distortion facors
RD_left = [k1_left, k2_left];
RD_right = [k1_right, k2_right];

% define intrinsic matrix from focal length and principal point
IM = [FL(1), 0, 0;
      0, FL(2), 0;
      PP(1), PP(2), 1];

% create a cameraParameters object for each camera
CP1 = cameraParameters('IntrinsicMatrix', IM, 'RadialDistortion', RD_left);
CP2 = cameraParameters('IntrinsicMatrix', IM, 'RadialDistortion', RD_right);

% manually define rotations around x, y, and z axes
Rx = [1,     0,       0    ;
       0,     cos(theta_x),    -sin(theta_x);
       0,  sin(theta_x),    cos(theta_x)   ];
   
Ry = [   cos(theta_y),    0, sin(theta_y);
          0,    1,    0    ;
       -sin(theta_y), 0,    cos(theta_y)   ];
   
Rz = [   cos(theta_z),  -sin(theta_z), 0;
       sin(theta_z),  cos(theta_z),     0;
          0,     0,     1];
      
R = Rz*Ry*Rx;

% create stereoParameters object and rectifc the images
stereoParams = stereoParameters(CP1,CP2,R,T);
[l_rect,r_rect] = rectifyStereoImages(l,r,stereoParams, 'OutputView', 'full');

%% Crop (by removing empty colums)
% first column with a value
crop = find(sum(rgb2gray(r_rect))~=0, 1, 'first');

% cropped to useful part
l_rect_crop = l_rect(:, 1:end-crop, :);
r_rect_crop = r_rect(:, crop+1:end, :);

%figure('name', "Cropped");
a = permute(cat(4, r_rect_crop, l_rect_crop), [4, 1, 2, 3]);
LFDispMousePan(a);

%% get sparse disparity map using feature correspondences
if strcmp(METHOD,'SPARSE')
    lg = im2single(rgb2gray(l_rect_crop));
    rg = im2single(rgb2gray(r_rect_crop));

    % enhance constrast with histogram equlization
    lg_eq = histeq(lg,255);
    rg_eq = histeq(rg, 255);

    % keypoint locations and descriptors

    %[Fl, Dl] = vl_sift(lg, 'PeakThresh', 0.001, 'EdgeThresh', 20);
    %[Fr, Dr] = vl_sift(rg, 'PeakThresh', 0.001, 'EdgeThresh', 20);


    % alternatively with dense SIFT at a single scale on a fixed grid
    % this is MUCH faster (speedup: 30-70x per descriptor)

    binSize = 7; % size of SIFT-descriptor bins
    step = [3, 10]; % extract a descriptor each STEP picels
    magnif = 3 ; % maginification factor to get smoothing right

    tic
    ls = vl_imsmooth(lg_eq, sqrt((binSize/magnif)^2 - .25)) ;
    rs = vl_imsmooth(rg_eq, sqrt((binSize/magnif)^2 - .25)) ;

    % perform dense descriptor extraction
    [Fl_dense, Dl_dense] = vl_dsift(ls, 'size', binSize, 'Step', step, 'Fast') ;
    [Fr_dense, Dr_dense] = vl_dsift(rs, 'size', binSize, 'Step', step, 'Fast') ;

    %% match the descriptors
    % rearange so that x runs first
    [~, perm] = sort(Fl_dense(2, :), 'ascend');
    Fl_dense = Fl_dense(:, perm);
    Fr_dense = Fr_dense(:, perm);
    Dl_dense = Dl_dense(:, perm);
    Dr_dense = Dr_dense(:, perm);

    % get number of features in each direction (y runs first)
    num_x = size(unique(Fl_dense(1, :)), 2);
    num_y = size(unique(Fl_dense(2, :)), 2);
    num_feats = size(Fl_dense, 2);

    % divide problem into horizontal subproblems and match
    thresh = 1.5; % default is 1.5. Greater -- less permissive

    matches = [];
    for i = 1:num_y
        indices = (i-1)*num_x+1:i*num_x;
        [matches_i, ~] = vl_ubcmatch(Dl_dense(:, indices), Dr_dense(:, indices), thresh);
        matches = [matches, matches_i+(i-1)*num_x];

    end

    %% Get Disparity from matches
    disp_range = [-75, 50];

    num_matches = size(matches, 2);
    matches_post = matches;
    disp = nan(1, num_matches);
    for i = 1:num_matches
        pos1 = Fl_dense(1:2, matches(1, i));
        pos2 = Fr_dense(1:2, matches(2, i));

        dx = pos1(1) - pos2(1);
        dy = pos1(2) - pos2(2);

        if dx < disp_range(1) || dx > disp_range(2)
            matches_post(:, i) = nan;
        else
            disp(i) = dx;
        end
    end

    % remove invalid matches from matches and disparity
    matches_post = matches_post(:, find(~isnan(matches_post(1, :))));
    disp = disp(:, find(~isnan(disp(1, :))));
    toc

    %figure();
    %histogram(disparity);

    %% Visualize GOOD matches as points on left image

    JL = im2uint8(l_rect_crop);

    xa = Fl_dense(1,matches_post(1,:)) ;
    ya = Fl_dense(2,matches_post(1,:)) ;

    pos = [xa', ya'];
    pos(:, 3) = 1;

    col_range = [-75, 50];
    col = (disp-col_range(1))/(col_range(2) - col_range(1));
    Color = uint8([1-col', 1-2*abs(col'-0.5), col']*255);


    J_circles = insertShape(JL,'circle',pos,'LineWidth',3, 'Color', Color);
    figure()
    imshow(J_circles);


    %% calculate REAL depth from sparse disparity
    Baseline = norm(T);
    real_disparity = crop + disp;

    z = FocalLength * Baseline ./ real_disparity; % in [mm]

    %% create x, y, z pointcloud (projection - orthographic)

    % get pixel positions of each point in z (use left image)
    x_pixel = Fl_dense(1,matches_post(1,:)) - size(JL,2)/2;
    y_pixel = Fl_dense(2,matches_post(1,:)) - size(JL,1)/2;

    x = z .* x_pixel ./ FL(1);
    y = z .* y_pixel ./ FL(2);

    ptCloud = pointCloud([x', y', z']);
    
    
%% Do the DENSE depthmap calculation
elseif strcmp(METHOD, 'DENSE')
    
    % define Range in which we look for disparity. This takes some tuning,
    % but in the future the process should be automeated (or we should make
    % an educated guess)
    disparityRange = round([-75, 50]/16)*16;
    
    lg = im2single(rgb2gray(l_rect_crop));
    rg = im2single(rgb2gray(r_rect_crop));

    % enhance constrast with histogram equlization
    lg_eq = histeq(lg,255);
    rg_eq = histeq(rg, 255);
    
    % get disparity Map
    disp = disparity(lg_eq, rg_eq, ...
                        'BlockSize',            15, ...
                        'ContrastThreshold',    0.5, ...
                        'UniquenessThreshold',  15, ...
                        'DisparityRange',       disparityRange);
                    
    % remove the edge values of the disparity range
    disp(disp>=disparityRange(2)-1) = -realmax('single');
    disp(disp<=disparityRange(1)) = -realmax('single');
    
    % Use some median filtering on the depthmap to smoothen more
    disp_filt = medfilt2(disp, [1, 1]);
    
    % remove tiny regions surrounded by much different values 
    % as they might be incorrect
    margin = 1;
    pixel_thresh = 1000;
    levels = disparityRange(1):0.5:disparityRange(2);
    complete_map = zeros(size(disp_filt));
    for i = 1:length(levels)
        l = levels(i);
        % get binary map
        bin_map = (disp_filt >=l-margin) & (disp_filt <= l+margin);
        CC = bwconncomp(bin_map);
        % remove components that are too small (< pixel_thresh)
        numPixels = cellfun(@numel,CC.PixelIdxList);
        idx = find(numPixels < pixel_thresh);
        for j = 1:length(idx)
            bin_map(CC.PixelIdxList{idx(j)}) = 0;
        end
        % add valid regions to overall map using logical OR
        complete_map = complete_map | bin_map;
    end
    
    % apply mask to disparity map
    complete_map = single(complete_map);
    complete_map(~complete_map) = nan;
    disp = complete_map.*disp_filt;
                    
    % display disparity map
    figure();
    imshow(disp_filt,disparityRange);
    title('Disparity Map from crop');
    colormap(gca,jet) 
    colorbar

    %% calculate REAL depth from sparse disparity
    Baseline = norm(T);    
    disp(disp<min(disparityRange)) = nan;
  
    real_disparity = crop + disp;

    z = FocalLength * Baseline ./ real_disparity; % in [mm]
    
    %% downsample disparity map to match sampling density of 3D model
    % 3D model x/y sampling is 0.195mm
    % here, one pixel corresponds approximately to avg(depth) / focallength
    z_mean = mean(z(find(~isnan(z))));
    downscale =  round(0.195 / (z_mean / FL(1)) / (2*sqrt(2))); % in [mm]
    z_down = z(1:downscale:end, 1:downscale:end);
    l_rect_crop_down = l_rect_crop(1:downscale:end, 1:downscale:end, :);

    %% create x, y, z pointcloud (projection - orthographic)

    % get pixel positions of each point in z (use left image)
    x_ax_pixel = (1:size(z, 2)) - size(z, 2)/2;
    y_ax_pixel = (1:size(z, 1)) - size(z, 1)/2;
    [x_pixel, y_pixel] = meshgrid(x_ax_pixel, y_ax_pixel);

    x = reshape(z .* x_pixel ./ FL(1), 1, []);
    y = reshape(z .* y_pixel ./ FL(2), 1, []);
    
        
    x_pixel = reshape(x_pixel, 1, []);
    y_pixel = reshape(y_pixel, 1, []);
    z = reshape(z, 1, []);
    
    colors = reshape(im2double(l_rect_crop), [], 3);

    ptCloud = pointCloud([x', y', z'], 'Color', colors);
    
    % get pixel positions of each point in z (use left image) for the
    % downsampeled depthmap
    xd_ax_pixel = ((1:size(z_down, 2)) - size(z_down, 2)/2) * downscale;
    yd_ax_pixel = ((1:size(z_down, 1)) - size(z_down, 1)/2) * downscale;
    [xd_pixel, yd_pixel] = meshgrid(xd_ax_pixel, yd_ax_pixel);

    xd = reshape(z_down .* xd_pixel ./ FL(1), 1, []);
    yd = reshape(z_down .* yd_pixel ./ FL(2), 1, []);
    
        
    xd_pixel = reshape(xd_pixel, 1, []);
    yd_pixel = reshape(yd_pixel, 1, []);
    zd = reshape(z_down, 1, []);
    
    colors = reshape(im2double(l_rect_crop_down), [], 3);

    ptCloud_down = pointCloud([xd', yd', zd'], 'Color', colors);
end
%% Visualization
if strcmp(METHOD, 'DENSE')
    figure();
    pcshow(ptCloud_down, 'MarkerSize', 60);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end

figure();
pcshow(ptCloud, 'MarkerSize', 12);
xlabel('X');
ylabel('Y');
zlabel('Z');