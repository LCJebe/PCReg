clear all
close all

%% define method for sampling points
% 'GRID' for regular grid
% 'RANDOM' for random points in point cloud
% 'ALL' for all points in point cloud
method = 'GRID';

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'Surface_DS_aligned.pcd'));
pcModel = pcread(strcat(path, 'GoodCrop.pcd'));
pcRand = pcread(strcat(path, 'RandCrop.pcd'));

%% define locations of spheres for local descriptors
if strcmp(method, 'GRID')
    XLimits = [min(pcSurface.XLimits(1), pcModel.XLimits(1)), ...
                max(pcSurface.XLimits(2), pcModel.XLimits(2))];
    YLimits = [min(pcSurface.YLimits(1), pcModel.YLimits(1)), ...
                max(pcSurface.YLimits(2), pcModel.YLimits(2))];
    ZLimits = [min(pcSurface.ZLimits(1), pcModel.ZLimits(1)), ...
                max(pcSurface.ZLimits(2), pcModel.ZLimits(2))];

    % create meshgrid of points
    d = 0.5;
    x = XLimits(1):d:XLimits(2);
    y = YLimits(1):d:YLimits(2);
    z = ZLimits(1):d:ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
end

%% get local descriptors for each model

% specify minimum number of points that has to be in sphere
min_pts = 30;
% specify radius of sphere
R = 1;

pts = pcModel.Location;
descModel = getLocalDescriptors(pts, sample_pts, min_pts, R);

pts = pcRand.Location;
descRand = getLocalDescriptors(pts, sample_pts, min_pts, R);

pts = pcSurface.Location;
descSurface = getLocalDescriptors(pts, sample_pts, min_pts, R);

%% function to calculate a descriptor for each point
function desc = getLocalDescriptors(pts, sample_pts, min_pts, R)
    % pts: points in pointcloud
    % sample_pts: points to calculate descriptors at
    % min_points: minimum number of points in sphere
    % R: Radius of sphere

    desc = [];
    tic
    for i = 1:size(sample_pts, 1)
        c = sample_pts(i, :);
        pts_local = getLocalPoints(pts, R, c, min_pts);

        if ~ isempty(pts_local) 
            num_points = size(pts_local, 1);

            % get first order moments [X, Y Z]
            M1_L2 = mean(pts_local, 1);
            M1_L1 = median(pts_local, 1);

            % get second order moments [XX, YY, ZZ, XY, XZ, YZ];
            M2 = (pts_local' * pts_local) / num_points;
            M2_L2 = [M2(1,1), M2(2,2), M2(3,3), M2(1,2), M2(1,3), M2(2, 3)];

            % get third order moments
            pts2 = pts_local.^2;
            % [XXX, YYY, ZZZ, XXY, XXZ, YYX, YYZ, ZZX, ZZY]
            M3 = (pts_local'.^2 * pts_local) / num_points;
            % [XYZ]
            mom3XYZ = mean(pts_local(:,1).*pts_local(:,2).*pts_local(:,3), 1);
            M3_L2 = [reshape(M3, 1, []), mom3XYZ];


            % descriptor structure is 
            % [X, Y, Z, X, Y, Z, XX, YY, ZZ, XY, XZ, YZ, XXX, YYY, ZZZ, XXY,
            % XXZ, YYX, YYZ, ZZX, ZZY, XYZ], 
            % where the first X, Y, Z are in L2 norm and the second X, Y, Z are
            % in L1 norm
            new_entry = [M1_L2, M1_L1, M2_L2, M3_L2];

            desc = [desc; new_entry];
        end
    end
    toc
end

%% Function: get points within sphere
% returns only the points from pts that are within a certain radius R
% of center point c if there are at least min_points in it
% RETURNS: Points RELATIVE to c. 
function pts_sphere = getLocalPoints(pts, R, c, min_points)

    % first: quickly crop cube around the center
    xLim = c(1) + [-R, R];
    yLim = c(2) + [-R, R];
    zLim = c(3) + [-R, R];
    mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
         & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
         & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
    mask = cat(2, mask, mask, mask);
    pts_cube = reshape(pts(mask), [], 3);    
    
    if size(pts_cube, 1) < min_points
        pts_sphere = [];
    else
        % then: get distances for each remaining point and only keep the points
        % that are within the radius
        pts_rel = pts_cube - c;
        dists = vecnorm(pts_rel, 2, 2);
        mask = dists < R;
        mask = cat(2, mask, mask, mask);
        pts_sphere = reshape(pts_rel(mask), [], 3);    
    end
    if size(pts_sphere, 1) < min_points
        pts_sphere = [];
    end
end