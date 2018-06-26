clear all
close all

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'c_Surface_down_aligned.pcd'));
pcModel = pcread(strcat(path, 'c_ModelCrop.pcd'));
pcRand = pcread(strcat(path, 'c_RandCrop.pcd'));

%% center both pointclouds so that everything is in positive octant
% center
pcSurface = centerPointCloud(pcSurface, 'positive');
pcModel = centerPointCloud(pcModel, 'positive');
pcRand = centerPointCloud(pcRand, 'positive');

%% define parameters
% number of voxels without overlap
num_voxels_xyz = 4;

% 1 for no overlap of voxels, 2 for 50%, 3 for 67% etc.
overlap = 2;

% minimum number of points per voxel so that calculations are done
min_points = 15;

%% calculate voxel sizes from parameters
size_complete = [pcSurface.XLimits(2), pcSurface.YLimits(2), pcSurface.ZLimits(2)];
size_voxels_xyz = size_complete ./ num_voxels_xyz;
total_num_voxels = overlap .* (num_voxels_xyz-1) + 1;

%% define scalar metrics used to create historgrams (prefix: sc_)
sc_median = [];
sc_mean = [];
sc_stdL1 = [];
sc_stdL1_9010 = [];
sc_ecc1 = [];
sc_ecc2 = [];

%% iterate over all voxels 
tic
for ix = 0:total_num_voxels-1
    xLim = [ix/overlap(1), 1+ix/overlap] * size_voxels_xyz(1);
    cx = xLim(1) + 0.5*size_voxels_xyz(1);
    for iy = 0:total_num_voxels-1
        yLim = [iy/overlap, 1+iy/overlap] * size_voxels_xyz(2);
        cy = yLim(1) + 0.5*size_voxels_xyz(2);
        for iz = 0:total_num_voxels-1
            zLim = [iz/overlap, 1+iz/overlap] * size_voxels_xyz(3);
            cz = zLim(1) + 0.5*size_voxels_xyz(3);
            
            % now we have the limits and the center of each voxel. Select
            % the points for each voxel
            pts = pcSurface.Location;
            mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
                 & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
                 & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
            mask = cat(2, mask, mask, mask);
            pts_vox = reshape(pts(mask), [], 3);    
                    
            % many voxels might be empty. check that at least minimum
            % number of points is included
            num_points = size(pts_vox, 1);
            if num_points > min_points
                pcVox = pointCloud(pts_vox);
                
                % center and align pointcloud with PCA
                pcRot = pcaCenterRotate(pcVox);
                pts_rot = pcRot.Location;
                
                % mean
                isc_mean = mean(pts_rot, 1);
                sc_mean = [sc_mean; isc_mean];
                
                % median
                isc_median = median(pts_rot, 1);
                sc_median = [sc_median; isc_median];
                
                % stdL1
                median_dists = pts_rot - isc_median;
                isc_stdL1 = median(median_dists, 1);
                sc_stdL1 = [sc_stdL1; isc_stdL1];
                
            end
        end
    end
end
toc

%% Function that centers a pointcloud
function pcCentered = centerPointCloud(PCin, method)
    % 'positive' sets the origin the the corner, so that all points are in
    % the positive octant
    if strcmp(method, 'positive')
        xc = PCin.XLimits(1);
        yc = PCin.YLimits(1);
        zc = PCin.ZLimits(1);

        A = eye(4);
        A(4, 1:3) = -[xc, yc, zc];
        tform = affine3d(A);

        pcCentered = pctransform(PCin,tform);
    % 'center' centers the pointcloud in the middle of extreme points
    elseif strcmp(method, 'center')
        xc = (PCin.XLimits(2) + PCin.XLimits(1))/2;
        yc = (PCin.YLimits(2) + PCin.YLimits(1))/2;
        zc = (PCin.ZLimits(2) + PCin.ZLimits(1))/2;

        A = eye(4);
        A(4, 1:3) = -[xc, yc, zc];
        tform = affine3d(A);

        pcCentered = pctransform(PCin,tform);
    end
end

%% Function that aligns pointcloud according to Principal components
function pcOut = pcaCenterRotate(PCin)
    % center pointcloud
    pcCent = centerPointCloud(PCin, 'center');
    
    % get PCA coefficients 
    points = pcCent.Location;
    coeff = pca(points);
    
    points_rot = points*coeff';
    pcOut = pointCloud(points_rot);
end