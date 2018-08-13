%% read in pointcloud
path = 'Data/PointClouds/';
pcName = 'Surface_DS2_aligned.pcd';
pcIn = pcread(strcat(path, pcName));

%% define translation offset and transform
t = [-5, -3, -1];
T = eye(4);
T(4, 1:3) = t;
tform = affine3d(T);

pcOut = pctransform(pcIn, tform);

%% Save aligned model to pcd file
save_file = strcat(path, pcName(1:end-4), '_offset.pcd');
pcwrite(pcOut, save_file);


%% same thing for a rotation
r = [0.5, 0.2, 0.8];
R = eul2rotm(r);

T = eye(4);
T(1:3, 1:3) = R;
tform = affine3d(T);

pcOut = pctransform(pcIn, tform);

%% Save aligned model to pcd file
save_file = strcat(path, pcName(1:end-4), '_rotated.pcd');
pcwrite(pcOut, save_file);