%% run getInliersRANSAC first to obtain a transformation T
% pcModel and pcSurface should also be in the workspace
% does not work yet!!!!!!!!!!!

%% apply estimated transformation to Surface (pts2)
% transform Surface from "raw" to "rand" location
load('TF_rand2raw.mat');
pts = quickTF(pcSurface.Location, invertTF(TF)); % invert so it's raw2rand

color = pcSurface.Color;
pts_tform = quickTF(pts, invertTF(T_final));

pcSurface_aligned = pointCloud(pts_tform, 'Color', color);

%% Save aligned model to pcd file
path = 'Data/PointClouds/';
%save_file = strcat(path, 'GoodCrop_autoalign.pcd');
%pcwrite(pcMode, save_file);
save_file = strcat(path, 'Surface_autoalign.pcd');
pcwrite(pcSurface_aligned, save_file);