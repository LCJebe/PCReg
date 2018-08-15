%% read in point cloud
pcModel = pcread('Data/PointClouds/GoodCropSmoothUp3.pcd');
pcRand = pcread('Data/PointClouds/RandCropSmoothUp3.pcd');

%% center / adjust
pcModel_centered = centerPointCloud(pcModel);
pcRand_centered = centerPointCloud(pcRand);

%% save pointcloud (careful: overwrite)
pcwrite(pcModel_centered, 'Data/PointClouds/GoodCropSmoothUp3.pcd');
pcwrite(pcRand_centered, 'Data/PointClouds/RandCropSmoothUp3.pcd');

%% helper function: center PC (0,0) to positive octant
function pc_out = centerPointCloud(pc)
    % translation of pointcloud
    T = [pc.XLimits(1), ...
        pc.YLimits(1), ...
        pc.ZLimits(1)];
    A = eye(4);
    A(4, 1:3) = -T;
    tform = affine3d(A);
    pc_out = pctransform(pc, tform);
end