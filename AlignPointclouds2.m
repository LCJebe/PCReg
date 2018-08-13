clear all
%% READ DATA: model and surface point clouds from saved file
path = 'Data/PointClouds/';
surface_name = 'Surface_DS2.pcd';
pcSurface = pcread(strcat(path, surface_name));
pcModel = pcread(strcat(path, 'GoodCrop.pcd'));

pcSurface = centerPointCloud(pcSurface);
pcModel = centerPointCloud(pcModel);

%% ALIGN: change coordinate system of surface

% specify rotation (in DEG) and translation (in mm)
x_off = 0.8; 
y_off = -0.8;
z_off = 1.3; 
Rx = 70; % 70
Ry = 188; % 188
Rz = 0; % 0

% build transformation matrix 
R = eul2rotm([Rx, Ry, Rz]/180*pi, 'XYZ');
T = [x_off, y_off, z_off];

A = eye(4);
A(1:3, 1:3) = R;
A(4, 1:3) = T;

% create tform object
tform = affine3d(A);

% transform pointcloud
pcAligned = pctransform(pcSurface,tform);

%% SAVE centered model 
save_file = strcat(path, 'GoodCrop_centered.pcd');
pcwrite(pcModel, save_file);

%% SAVE resulting pointcloud surface as pcd
save_file = strcat(path, surface_name(1:end-4), '_aligned_centered.pcd');
pcwrite(pcAligned, save_file);

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