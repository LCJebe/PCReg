%% run getInliersRANSAC first to load necessary variables into workspace

path = 'Data/PointClouds/';
%% create one pointcloud that contains matched point pairs post-alignment

color1 = repmat([255, 0, 0], size(inlierPtIdx, 1), 1); % RED: Model
color2 = repmat([0, 255, 0], size(inlierPtIdx, 1), 1); % GREEN: Surface
color = [color1; color2];

pts = single([pts1_aligned(inlierPtIdx, :); pts2(inlierPtIdx, :)]);
pcInliers = pointCloud(pts, 'Color', color);


%% save Inliers
save_file = strcat(path, 'inliers_aligned_DELETE.pcd');
pcwrite(pcInliers, save_file);

%% DEBUG: Save all matches
pts = single([pts1_aligned; pts2]);

color1 = repmat([255, 0, 0], size(pts, 1)/2, 1); % RED: Model
color2 = repmat([0, 255, 0], size(pts, 1)/2, 1); % GREEN: Surface
color = [color1; color2];

pcInliers = pointCloud(pts, 'Color', color);

save_file = strcat(path, 'matches_aligned_DELETE.pcd');
pcwrite(pcInliers, save_file);

%% DEBUG: Save matches after manual alignment with align2

pts = single([pts1; pts2_postalign]);

color1 = repmat([255, 0, 0], size(pts, 1)/2, 1); % RED: Model
color2 = repmat([0, 255, 0], size(pts, 1)/2, 1); % GREEN: Surface
color = [color1; color2];

pcInliers = pointCloud(pts, 'Color', color);

save_file = strcat(path, 'matches_postalign_DELETE.pcd');
pcwrite(pcInliers, save_file);