%% read in pointcloud
path = 'Data/PointClouds/';
pcName = 'Surface_DS_aligned.pcd';
pcIn = pcread(strcat(path, pcName));

%% define offset and transform
t = [-5, -3, -1];
T = eye(4);
T(4, 1:3) = t;
tform = affine3d(T);

pcOut = pctransform(pcIn, tform);

%% Save aligned model to pcd file
save_file = strcat(path, pcName(1:end-4), '_offset.pcd');
pcwrite(pcOut, save_file);