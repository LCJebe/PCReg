%% run getInliersRANSAC first to obtain a transformation T
% pcModel and pcSurface should also be in the workspace

%% apply estimated transformation to Model (pts1)
tform = affine3d(T);

pcModel_aligend = pctransform(pcModel, tform);

%% Save aligned model to pcd file
path = 'Data/PointClouds/';
save_file = strcat(path, 'GoodCrop_autoalign.pcd');
pcwrite(pcModel_aligend, save_file);
save_file = strcat(path, 'Surface_autoalign.pcd');
pcwrite(pcSurface, save_file);