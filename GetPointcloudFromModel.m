clear all
close all

% 0: cropped to full model
% 1: manual box select
% 2: predefined settings for good crop
% 3: predefined settings for bad "random" crop
% 4: predefined settings for 2nd bad "random" crop
% 5: predefined settings for 3rd bad "random" crop
% 6: predefined settings for half-good, half-bad crop
% 7: no crop
CROP = 0;

enable_3D_edges = true;

%% Read DICOM images
% 512 x 512
%dicom_path = 'Data/DICOM/Felsenbein/DICOM512/';
% 1024 x 1024
dicom_path = 'Data/DICOM/Felsenbein/DICOM1024/';
% Cadaver Scan
%dicom_path = 'Data/DICOM/Cadaver_Postop/';

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

%% Rescale to unit8 using range [0, 4095], remove singleton channel dimension
% order of axes: y, x, z
minV = min(V(:));
maxV = max(V(:));
rangeV = maxV - minV + 1;
scale_factor = double(rangeV) / 4096;
Vs = V / scale_factor;
Vs = Vs - min(Vs(:));
VV = uint8(double(squeeze(Vs))*255/4096);
[m, n, d] = size(VV);

if CROP == 0
    ymin = 444;
    ymax = 729;
    xmin = 313;
    xmax = 832;
    zmin = 12;
    zmax = 302;
elseif CROP == 1
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
elseif CROP == 4 % incorrect area sampeled
    ymin = 460;
    ymax = 530;
    xmin = 440;
    xmax = 540;
    zmin = 60;
    zmax = 100;
elseif CROP == 5 % incorrect area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 440;
    xmax = 540;
    zmin = 20;
    zmax = 60;
elseif CROP == 6 % challenge: half-good, half-bad area sampeled
    ymin = 530;
    ymax = 600;
    xmin = 490;
    xmax = 590;
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
figure('name', 'Cropped Volume');
LFDispMousePan(permute(Vc, [3, 1, 2]), 2);


%% Binarize Image and get rid of unneccesary parts
[mc, nc, dc] = size(Vc);

% threshold for density: 0.25 (==1000) should be air...
thresh = 0.2; % as percentage between 0.2 (air?) and 4000 (metal)

%Vb: binarized Volume
Vb = zeros(size(Vc));
for iz = 1:dc
    % binarize slice with given threshold value
    slice = squeeze(Vc(:, :, iz));
    slice_bin = imbinarize(slice, thresh);
    Vb(:, :, iz) = slice_bin;
end

% show binary volume
figure('name', 'Binary Volume');
LFDispMousePan(permute(Vb, [3, 1, 2]), 2);


%% if desired, fill holes to reduce surface area
FILL_HOLES = true;
min_pixels = round(mc*nc*dc/5); % ideally, there is only one single air region!

if FILL_HOLES
    r = 5;   
    % expand volume with zeros to not remove holes that appear at edges
    Vb_expanded = padarray(Vb,[r, r, r],0,'both');
    % remove small regions with 6-connectivity for speed
    Vb_expanded = ~bwareaopen(~Vb_expanded, min_pixels, 6);
    Vb = double(Vb_expanded(r+1:end-r, r+1:end-r, r+1:end-r));
end

% show binary volume after hole removal
figure('name', 'Binary Volume, holes filled');
LFDispMousePan(permute(Vb, [3, 1, 2]), 2);

%% Edge detection 3D

% create custom colormap
span = 256 - round(thresh*255);
cmap = jet(span);
custom_cmap = zeros(256, 3);
custom_cmap(end-span+1:end, :) = cmap;


if enable_3D_edges
    
    % detect all edges in 3D    
%     kernel_6conn = ones(3, 3, 3);
%     kernel_6conn(1, 2, 2) = 21;
%     kernel_6conn(2, 1, 2) = 21;
%     kernel_6conn(2, 2, 1) = 21;
%     kernel_6conn(3, 2, 2) = 21;
%     kernel_6conn(2, 3, 2) = 21;
%     kernel_6conn(2, 2, 3) = 21;
%     kernel_6conn(2, 2, 2) = 147;    
% 
    kernel_6conn = zeros(3, 3, 3);
    kernel_6conn(1, 2, 2) = 1;
    kernel_6conn(2, 1, 2) = 1;
    kernel_6conn(2, 2, 1) = 1;
    kernel_6conn(3, 2, 2) = 1;
    kernel_6conn(2, 3, 2) = 1;
    kernel_6conn(2, 2, 3) = 1;
    kernel_6conn(2, 2, 2) = 7; 
    

    conMap6 = convn(Vb, kernel_6conn, 'valid');
    %conMap6 = padarray(conMap6, [1, 1, 1], 294, 'both'); 
    conMap6 = padarray(conMap6, [1, 1, 1], 13, 'both'); 
    Ve3D_6conn = Vb;
    %Ve3D_6conn(conMap6 >= 293) = 0;
    Ve3D_6conn(conMap6 == 13) = 0;

    
    % show edge Volume
    figure('name', '3D Detected Edges');
    LFDispMousePan(permute(Ve3D_6conn, [3, 1, 2, 4]), 2);
    
    % 3D Dilated Volume
    Vdil = imdilate(Vc, ones(3, 3, 3));
    % Density Colored Volume (slice by slice...)
    Vcol = zeros(mc, nc, dc, 3);
    for iz = 1:dc
        Vcol(:, :, iz, :) = ind2rgb(squeeze(Vdil(:, :, iz)), custom_cmap);
    end
    % Density colored edges
    Ve3D_6conn_col = uint8(double(Ve3D_6conn).*Vcol*255);
    
    % show edge Volume
    figure('name', 'Density Colored Surface Edges');
    LFDispMousePan(permute(Ve3D_6conn_col, [3, 1, 2, 4]), 2);
else
%% Edge detection 2D (boundaries of binary image)

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
    figure('name', 'Density Colored Surface Edges');
    LFDispMousePan(permute(Vcolor_edges, [3, 1, 2, 4]), 2);
end

%% create pointcloud from detected edges

% define arbitrary 0,0,0 point at position 1,1,1 and use dy, dx, dz to
% convert to real world coordinates in [mm]


if enable_3D_edges
    [yp, xp, zp] = ind2sub(size(Ve3D_6conn), find(Ve3D_6conn));
    Ve_expanded = cat(4, Ve3D_6conn, Ve3D_6conn, Ve3D_6conn); 
    colors = reshape(Ve3D_6conn_col(logical(Ve_expanded)), [], 3);
else
    [yp, xp, zp] = ind2sub(size(Ve), find(Ve));
    Ve_expanded = cat(4, Ve, Ve, Ve); 
    colors = reshape(Vcolor_edges(logical(Ve_expanded)), [], 3);
end

xyz_points = single([xp, yp, zp] .* [dx, dy, dz]);

ptCloud = pointCloud(xyz_points, 'Color', colors);  % with real [mm] scaling
%ptCloud = pointCloud(single([-xp, zp, yp])); % equidistant in samples
num_points = size(xyz_points, 1);

%% Visualization
figure('name', 'Point Cloud');
pcshow(ptCloud, 'MarkerSize', 30, 'VerticalAxis', 'Y');
xlabel('X');
ylabel('Y');
zlabel('Z');
