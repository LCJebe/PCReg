%% run getInliersRANSAC first to obtain a transformation T
% pcModel and pcSurface should also be in the workspace
% does not work yet!!!!!!!!!!!

%% apply estimated transformation to Surface (pts2)

% load and center the correct surface
% raw = 'SurfaceNew_DS3.pcd'
% rand = 'SurfaceNew_DS3.pcd'
% rand2 = 'SurfaceNew_DS3_rand2.pcd'
pcSurface = pcread('Data/PointClouds/SurfaceNew_DS3.pcd');
pcSurface = centerPointCloud(pcSurface);


USE_RAND2RAW = false; % true if the workspace ends with _rand.mat
if USE_RAND2RAW
    % transform Surface from "raw" to "rand" location
    load('TF_rand2raw.mat');
    pts = quickTF(pcSurface.Location, invertTF(TF)); % invert so it's raw2rand
else
    pts = pcSurface.Location;
end

color = pcSurface.Color;
pts_tform = quickTF(pts, invertTF(T_final));

pcSurface_aligned = pointCloud(pts_tform, 'Color', color);

%% Save aligned model to pcd file
path = 'Data/PointClouds/';
%save_file = strcat(path, 'GoodCrop_autoalign.pcd');
%pcwrite(pcMode, save_file);
save_file = strcat(path, 'Surface_autoalign.pcd');
pcwrite(pcSurface_aligned, save_file);


%% helper function: center PC (0,0) in middle of crop
function pc_out = centerPointCloud(pc)
    % translation of pointcloud
    T = [(pc.XLimits(1) + pc.XLimits(2))/2, ...
        (pc.YLimits(1) + pc.YLimits(2))/2, ...
        (pc.ZLimits(1) + pc.ZLimits(2))/2];
    A = eye(4);
    A(4, 1:3) = -T;
    tform = affine3d(A);
    pc_out = pctransform(pc, tform);
end