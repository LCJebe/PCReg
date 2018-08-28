clear all
%% READ DATA: model and surface point clouds from saved file
path = 'Data/PointClouds/';
surface_name = 'SurfaceNew_DS3.pcd';
reference_name = 'Surface_DS2_alignedM.pcd';
pcSurface = pcread(strcat(path, surface_name));
pcReference = pcread(strcat(path, reference_name));

%% ALIGN: change coordinate system of surface
% center first 
pcSurface = centerPointCloud(pcSurface);

% specify rotation (in DEG) and translation (in mm)
x_off = 35.2; % 11.7
y_off = 22.8; % 6.1
z_off = 24.5; % 8.3
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

% transform pointcloud (rotate) 
pcAligned = pctransform(pcSurface,tform);

%% now that the point clouds are rotated nicely, get limits and align
AUTOTRANSLATE = false;
if AUTOTRANSLATE
    T = [pcReference.XLimits(1) - pcAligned.XLimits(1), ...
         pcReference.YLimits(1) - pcAligned.YLimits(1), ...
         pcReference.ZLimits(1) - pcAligned.ZLimits(1)];

    A = eye(4);
    A(4, 1:3) = T;
    tform = affine3d(A);

    % transform pointcloud (translate)
    pcAligned = pctransform(pcAligned, tform);
end

%% SAVE resulting pointcloud surface as pcd
save_file = strcat(path, surface_name(1:end-4), '_alignedM.pcd');
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