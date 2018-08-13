%% run getInliersRANSAC first to obtain a transformation T
% pcModel and pcSurface should also be in the workspace

%% apply estimated transformation to Model (pts1)

pts = pcModel.Location;
color = pcModel.Color;
pts_tform = [pts, ones(size(pts, 1), 1)] * T;

pcModel_aligend = pointCloud(pts_tform(:, 1:3), 'Color', color);

%% Save aligned model to pcd file
path = 'Data/PointClouds/';
save_file = strcat(path, 'GoodCrop_autoalign.pcd');
pcwrite(pcModel_aligend, save_file);
save_file = strcat(path, 'Surface_autoalign.pcd');
pcwrite(pcSurface, save_file);