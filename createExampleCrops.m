%% load full model
pc = pcread('Data/PointClouds/ModelSmoothUp3.pcd');

%% previous crop data (needs to be added first) this is the data
% on which the loaded model is based. 
ymin_pre = 444;
ymax_pre = 729;
xmin_pre = 313;
xmax_pre = 832;
zmin_pre = 12;
zmax_pre = 302;


%% select crop area
% 2: predefined settings for good crop
% 3: predefined settings for bad "random" crop
% 4: predefined settings for 2nd bad "random" crop
% 5: predefined settings for 3rd bad "random" crop
% 6: predefined settings for half-good, half-bad crop
CROP = 3;

% define name for saving
savename = 'RandCropSmoothUp3.pcd';

% specify crop box
if CROP == 2 % correct area centered
    ymin = 530;
    ymax = 600;
    xmin = 440;
    xmax = 540;
    zmin = 60;
    zmax = 100;
elseif CROP == 3 % incorrect area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 340;
    xmax = 440;
    zmin = 60;
    zmax = 100;
elseif CROP == 4 % incorrect area sampeled
    ymin = 460;
    ymax = 530;
    xmin = 440;
    xmax = 540;
    zmin = 60;
    zmax = 100;
elseif CROP == 5 % incorrect area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 440;
    xmax = 540;
    zmin = 20;
    zmax = 60;
elseif CROP == 6 % challenge: half-good, half-bad area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 490;
    xmax = 590;
    zmin = 60;
    zmax = 100;
end

% define world units and convert to world units
dx = 0.1953;
dy = 0.1953;
dz = 0.3400;

ymin = (ymin-ymin_pre)*dy;
ymax = (ymax-ymin_pre)*dy;
xmin = (xmin-xmin_pre)*dx;
xmax = (xmax-xmin_pre)*dx;
zmin = (zmin-zmin_pre)*dz;
zmax = (zmax-zmin_pre)*dz;

%% crop
pts = pc.Location;

mask = pts(:, 1) > xmin & pts(:, 1) < xmax ...
     & pts(:, 2) > ymin & pts(:, 2) < ymax ...
     & pts(:, 3) > zmin & pts(:, 3) < zmax;
mask = cat(2, mask, mask, mask);
pts_crop = reshape(pts(mask), [], 3);    

pc_crop = pointCloud(pts_crop);

%% optional: save pointCloud
pcwrite(pc_crop, strcat('Data/PointClouds/', savename));