clear all
close all

% 0: no crop
% 1: manual box select
% 2: predefined settings for good crop
% 3: predefined settings for bad "random" crop
CROP = 1;

%% Read DICOM images
% 512 x 512
%dicom_path = '/media/lars/Elements/ARRI/Felsenbein/Inn_Ear_512_i2_iDose_(1)_202/';
% 1024 x 1024
dicom_path = '/media/lars/Elements/ARRI/Felsenbein/Inn_Ear_1024_iDose_(4)_201/';

% read in data
%      V: 4D Volume as (rows, columns, channels, slices)
%      spatial: information about real world coordinates / units [mm]
%      dim: dimension with largest offset from previous slice (1, 2, 3)
[V,spatial,dim] = dicomreadVolume(dicom_path);

% get spatial offset in x, y, z direction, assuming they are constant
% across the whole model, all in [mm]
dx = spatial.PixelSpacings(1,1);
dy = spatial.PixelSpacings(1,2);
dz = spatial.PatientPositions(2, 3) - spatial.PatientPositions(1, 3);

%% Rescale to unit8 using range [0, 3000], remove channel dimension
% order of axes: y, x, z
VV = uint8(double(squeeze(max(min(V, 3000), 0)))*255/3000);
[m, n, d] = size(VV);

if CROP == 1
    % display center in y, x, z view and let user select boundaries
    figy = figure();
    imshow(squeeze(VV(uint16(m/2), :, :)));
    recty = getrect(figy);
    figx = figure();
    imshow(squeeze(VV(:, uint16(n/2), :)));
    rectx = getrect(figx);
    figz = figure();
    imshow(squeeze(VV(:, :, uint16(d/2))));
    rectz = getrect(figz);
    close all

    % get largest selection for x, y, z and crop accordingly

    xmin = uint16(max(min(rectz(1), recty(2)), 1));
    xmax = uint16(min(max(rectz(1)+rectz(3), recty(2)+recty(4)), n));
    ymin = uint16(max(min(rectz(2), rectx(2)), 1));
    ymax = uint16(min(max(rectz(2)+rectz(4), rectx(2)+rectx(4)), m));
    zmin = uint16(max(min(rectx(1), recty(1)), 1));
    zmax = uint16(min(max(rectx(1)+rectx(3), recty(1)+recty(3)), d));
    
elseif CROP == 2 % correct area centered
    ymin = 530;
    ymax = 600;
    xmin = 440;
    xmax = 540;
    zmin = 60;
    zmax = 100;
elseif CROP == 3 % incorrect area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 340;
    xmax = 440;
    zmin = 60;
    zmax = 100;
else
    xmin = 1;
    xmax = n;
    ymin = 1;
    ymax = m;
    zmin = 1;
    zmax = d;
end

% order: (y, x, z)
% Vc: cropped Volume
Vc = VV(ymin:ymax, xmin:xmax, zmin:zmax);

%show cropped volume
%figure();
%LFDispMousePan(permute(Vc, [3, 1, 2]));


%% Binarize Image and get rid of unneccesary parts
[mc, nc, dc] = size(Vc);

thresh = 0.3;


%Vb: binarized Volume
Vb = zeros(size(Vc));
for iz = 1:dc
    % binarize slice with given threshold value
    slice = squeeze(Vc(:, :, iz));
    slice_bin = imbinarize(slice, thresh);

    Vb(:, :, iz) = slice_bin;
end

% show binary volume
figure();
LFDispMousePan(permute(Vb, [3, 1, 2]));


%% if desired, fill holes to reduce surface area!!
FILL_HOLES = true;
min_pixels = 1e6; % ideally, there is only one single air region!

if FILL_HOLES
    r = 5;   
    % expand volume with zeros to not remove holes that appear at edges
    Vb_expanded = padarray(Vb,[r, r, r],0,'both');
    % remove small regions with 6-connectivity for speed
    Vb_expanded = ~bwareaopen(~Vb_expanded, min_pixels, 6);
    Vb = double(Vb_expanded(r+1:end-r, r+1:end-r, r+1:end-r));
end

% show binary volume after hole removal
figure();
LFDispMousePan(permute(Vb, [3, 1, 2]));


%% Edge detection

% create custom colormap
span = 256 - round(thresh*255);
cmap = jet(span);
custom_cmap = zeros(256, 3);
custom_cmap(end-span+1:end, :) = cmap;

%Ve: Edge Volume
Ve = zeros(size(Vb), 'uint8');
Vcolor_edges = zeros(mc, nc, dc, 3, 'uint8');
for iz = 1:dc
    slice_dil = imdilate(squeeze(Vc(:, :, iz)), ones(7));
    slice_dil_falsecolor = ind2rgb(slice_dil, custom_cmap);
    slice_bin = squeeze(Vb(:, :, iz));
    % detect edges on a slice using canny
    B = bwboundaries(slice_bin);
    edges_img = zeros(size(slice_bin), 'uint8');
    for k = 1:length(B)
        boundary = B{k};
        edges_img(sub2ind(size(slice_bin), boundary(:, 1), boundary(:, 2))) = 255;
        % dont consider image boundaries as edges!
        edges_img(1, :) = 0;
        edges_img(:, 1) = 0;
        edges_img(end, :) = 0;
        edges_img(:, end) = 0;
    end
    Ve(:, :, iz) = edges_img;
    Vcolor_edges(:,:,iz,:) = uint8(round(im2double(edges_img).*slice_dil_falsecolor*255));
end

% show edge Volume
%figure();
%LFDispMousePan(permute(Vcolor_edges, [3, 1, 2, 4]), 2);

%% create pointcloud from detected edges

% define arbitrary 0,0,0 point at position 1,1,1 and use dy, dx, dz to
% convert to real world coordinates in [mm]

[yp, xp, zp] = ind2sub(size(Ve), find(Ve));
Ve_expanded = cat(4, Ve, Ve, Ve); 
colors = reshape(Vcolor_edges(logical(Ve_expanded)), [], 3);
xyz_points = single([xp, yp, zp] .* [dx, dy, dz]);

ptCloud = pointCloud(xyz_points, 'Color', colors);  % with real [mm] scaling
%ptCloud = pointCloud(single([-xp, zp, yp])); % equidistant in samples
num_points = size(xyz_points, 1);

%% Visualization
figure();
pcshow(ptCloud, 'MarkerSize', 30, 'VerticalAxis', 'Y');
xlabel('X');
ylabel('Y');
zlabel('Z');
