% define pointcloud
num_pts = 1e5;
rangeX = 10;
rangeY = 10;
rangeZ= 10;

% create random pointcloud
sample_pts = rand(num_pts, 3);

% scale numbers so that they fit into the correct range
sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];

pc = pointCloud(sample_pts);
pcshow(pc);

% define center and radius and get local points
c = [3, 3, 3];
R = 2.5;

ptsLocal = getLocalPoints(sample_pts, R, c, 100);
figure()
pcshow(pointCloud(ptsLocal));