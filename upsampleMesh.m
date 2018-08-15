%% load the mesh into workspace
[vertices, faces] = read_ply('Data/Mesh/ModelSmooth.ply');

%% upsample points by creating one additional point per triangle face

% each pt vector (n x 1) specifies one corner of a triangle for all triangles
pt1 = faces(:, 1);
pt2 = faces(:, 2);
pt3 = faces(:, 3);

% each pos vector (n x 3) holds the coordinates of the pt vector points
pos1 = vertices(pt1, :);
pos2 = vertices(pt2, :);
pos3 = vertices(pt3, :);

% to create a new point on the triangle surface, average the three pos
% vectors
new_pts =  (pos1 + pos2 + pos3) / 3;

%% save the old and new points as new pointcloud
pts = single([vertices; new_pts]);
clear new_pts pos1 pos2 pos3 pt1 pt2 pt3

pcwrite(pointCloud(pts), 'Data/PointClouds/ModelSmoothUp3.pcd');

