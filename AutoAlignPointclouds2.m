%% run getInliersRANSAC first to obtain a transformation T
% pcModel and pcSurface should also be in the workspace
% does not work yet!!!!!!!!!!!

%% apply estimated transformation to Surface (pts2)

pts = pcSurface.Location;
color = pcSurface.Color;
pts_tform = [pts, ones(size(pts, 1), 1)]*T;

pcSurface_aligned = pointCloud(pts_tform(:, 1:3), 'Color', color);

%% Save aligned model to pcd file
path = 'Data/PointClouds/';
%save_file = strcat(path, 'GoodCrop_autoalign.pcd');
%pcwrite(pcMode, save_file);
save_file = strcat(path, 'Surface_autoalign.pcd');
pcwrite(pcSurface_aligned, save_file);