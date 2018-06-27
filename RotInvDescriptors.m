clear all
close all

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'c_Surface_down_aligned.pcd'));
pcModel = pcread(strcat(path, 'c_ModelCrop.pcd'));
pcRand = pcread(strcat(path, 'c_RandCrop.pcd'));



%% define parameters
% voxel size (cubic voxels) in [mm]
params.voxel_size = 3;

% 1 for no overlap of voxels, 2 for 50%, 3 for 67% etc.
params.overlap = 2;

% minimum number of points per voxel so that calculations are done
params.min_points = 50;

%% get descriptor for each pointcloud
descSurface =  getDescriptor(pcSurface, params);
descModel =  getDescriptor(pcModel, params);
descRand =  getDescriptor(pcRand, params);

%% plot histograms
Nbins = 16;
hist_edges.mean = linspace(0.5, 1.5, Nbins+1);
hist_edges.median = linspace(0.5, 1.5, Nbins+1);
hist_edges.std = linspace(0, 0.4, Nbins+1);
hist_edges.stdL1 = linspace(0, 0.3, Nbins+1);
hist_edges.ecc1 = linspace(0, 1.5, Nbins+1);
hist_edges.ecc2 = linspace(0, 1.5, Nbins+1);

histsSurface = plotHistograms(descSurface, hist_edges, 'Surface');
histsModel = plotHistograms(descModel, hist_edges, 'Good Crop');
histsRand = plotHistograms(descRand, hist_edges, 'Random Crop');

%% get distance from surface to each histogram
[dGood, dBad, dGoodBad, dProof] = deal(zeros(1, 6));

for i = 1:6
    dGood(i) = histogram_intersection(histsSurface(i,:), histsModel(i,:));
    dBad(i) = histogram_intersection(histsSurface(i,:), histsRand(i,:));
    dGoodBad(i) = histogram_intersection(histsModel(i,:), histsRand(i,:));
    dProof = histogram_intersection(histsModel(i,:), histsModel(i,:));
end

fprintf('Distance to bad: %f\n', sum(dBad));
fprintf('Distance to good: %f\n', sum(dGood));
fprintf('Distance between crops: %f\n', sum(dGoodBad));
fprintf('Distance to self: %f\n', sum(dProof));

%% Function that returns a descriptor of size Nx6 for a given input volume
function desc = getDescriptor(Pcin, params)
    voxel_size = params.voxel_size;
    overlap = params.overlap;
    min_points = params.min_points;

    %% align center all pointclouds so that everything is in positive octant
    % this centers around the middle and aligns with PCA
    Pcin = pcaCenterRotate(Pcin);
    
    % this shift everything into the positive octant
    Pcin = centerPointCloud(Pcin, 'positive');
    
    %% calculate voxel sizes from parameters
    % only look at voxels that fit in there completely
    size_complete = [Pcin.XLimits(2), Pcin.YLimits(2), Pcin.ZLimits(2)];
    integer_num_voxels_xyz = floor(size_complete ./ voxel_size * overlap);
    total_num_voxels = integer_num_voxels_xyz-overlap + 1;

    %% define scalar metrics used to create historgrams (prefix: sc_)
    sc_median = [];
    sc_mean = [];
    sc_std = [];
    sc_stdL1 = [];
    sc_ecc1 = [];
    sc_ecc2 = [];

    %% iterate over all voxels 
    tic
    fprintf('Dividing pointcloud into %d x %d x %d voxel grid...\n', ...
            total_num_voxels(1), total_num_voxels(2), total_num_voxels(3));
    for ix = 0:total_num_voxels(1)-1
        xLim = [ix/overlap, 1+ix/overlap] * voxel_size;
        cx = xLim(1) + 0.5*voxel_size;
        for iy = 0:total_num_voxels(2)-1
            yLim = [iy/overlap, 1+iy/overlap] * voxel_size;
            cy = yLim(1) + 0.5*voxel_size;
            for iz = 0:total_num_voxels(3)-1
                zLim = [iz/overlap, 1+iz/overlap] * voxel_size;
                cz = zLim(1) + 0.5*voxel_size;

                % now we have the limits and the center of each voxel. Select
                % the points for each voxel
                pts = Pcin.Location;
                mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
                     & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
                     & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
                mask = cat(2, mask, mask, mask);
                pts_vox = reshape(pts(mask), [], 3);    

                % many voxels might be empty. check that at least minimum
                % number of points is included
                num_points = size(pts_vox, 1);
                if num_points > min_points

                    % mean distance from center
                    pts_rel = pts_vox - [cx, cy, cz];
                    pts_dists = vecnorm(pts_rel, 2, 2);
                    isc_mean = mean(pts_dists);
                    sc_mean = [sc_mean; isc_mean];

                    % median distance from center
                    isc_median = median(pts_dists);
                    sc_median = [sc_median; isc_median];

                    % std of distances from center
                    isc_std = std(pts_dists);
                    sc_std = [sc_std, isc_std];

                    % stdL1 of distances from center
                    median_dists = abs(pts_dists - isc_median);
                    isc_stdL1 = median(median_dists, 1);
                    sc_stdL1 = [sc_stdL1; isc_stdL1];

                    % for eccentricities, align each voxel via PCA
                    inVox = pointCloud(pts_vox);
                    rotVox = pcaCenterRotate(inVox);
                    % find difference between 90% and 10% L1 projections
                    dists_rot = abs(rotVox.Location);
                    i90 = ceil(0.9*num_points);
                    i10 = ceil(0.1*num_points);
                    % sort for each projection separately
                    x_sorted = sort(dists_rot(:, 1));
                    y_sorted = sort(dists_rot(:, 2));
                    z_sorted = sort(dists_rot(:, 3));
                    dx = x_sorted(i90) - x_sorted(i10);
                    dy = y_sorted(i90) - y_sorted(i10);
                    dz = z_sorted(i90) - z_sorted(i10);

                    isc_ecc1 = dy / dx;
                    isc_ecc2 = dz / dx;

                    sc_ecc1 = [sc_ecc1, isc_ecc1];
                    sc_ecc2 = [sc_ecc2, isc_ecc2];
                end
            end
        end
    end  
    fprintf('Found %d non-empty voxels!\n', size(sc_mean, 1));
    toc
    
    % normalize distances by the voxel size
    desc.mean = sc_mean / (0.5*voxel_size);
    desc.median = sc_median / (0.5*voxel_size);
    desc.std = sc_std / (0.5*voxel_size);
    desc.stdL1 = sc_stdL1 / (0.5*voxel_size);
    desc.ecc1 = sc_ecc1;
    desc.ecc2 = sc_ecc2;
end
    

%% Function that shows all 6 histograms in one plot
function histograms = plotHistograms(desc, hist_edges, win_title)
    figure('name', win_title);
    subplot(2,3,1);
    hists.mean = histogram(desc.mean, hist_edges.mean, 'Normalization', 'probability');
    title('Mean Dist to Center');


    subplot(2,3,2);
    hists.median = histogram(desc.median, hist_edges.median, 'Normalization', 'probability');
    title('Median Dist to Center');


    subplot(2,3,3);
    hists.std = histogram(desc.std, hist_edges.std, 'Normalization', 'probability');
    title('Std of Dist to Center');


    subplot(2,3,4);
    hists.stdL1 = histogram(desc.stdL1, hist_edges.stdL1, 'Normalization', 'probability');
    title('L1-Std of Dist to Center');


    subplot(2,3,5);
    hists.ecc1 = histogram(desc.ecc1, hist_edges.ecc1, 'Normalization', 'probability');
    title('Eccentricity 1');


    subplot(2,3,6);
    hists.ecc2 = histogram(desc.ecc2, hist_edges.ecc2, 'Normalization', 'probability');
    title('Ecentricity 2');
    
    % merge histograms into one matrix (each row is one histogram)
    histograms = zeros(6, length(hists.ecc2.BinCounts));
    histograms(1, :) = hists.mean.BinCounts / sum(hists.mean.BinCounts);
    histograms(2, :) = hists.median.BinCounts / sum(hists.median.BinCounts);
    histograms(3, :) = hists.std.BinCounts / sum(hists.std.BinCounts);
    histograms(4, :) = hists.stdL1.BinCounts / sum(hists.stdL1.BinCounts);
    histograms(5, :) = hists.ecc1.BinCounts / sum(hists.ecc1.BinCounts);
    histograms(6, :) = hists.ecc2.BinCounts / sum(hists.ecc2.BinCounts);
end

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
    
    points_rot = points*coeff; %% don't transpose the coefficients!!
    pcOut = pointCloud(points_rot);
end