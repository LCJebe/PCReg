%% load the ply file with vertices and faces
[pts_xyz, fcs] = read_ply('Data/Mesh/ModelSmooth.ply');

%% transform points into viewing coordinates
% Clipping plane [near,far] 69.7784, 307.341
% Focal point [x,y,z] 49.6397, 27.9948, 49.3149
% Position [x,y,z] 23.5838, -130.702, 117.717
% View up [x,y,z] 0.692896, -0.378331, -0.613808
% Camera view angle [degrees] 30
% Window size [x,y] 1855, 1056
% Window position [x,y] 65, 24

% Clipping Planes
cp = [69.7784, 307.341]; % near, far
%Focal Point
fp = [49.6397, 27.9948, 49.3149]; % x, y, 
% View Reperence Point (Camera Position)
VRP = [23.5838, -130.702, 117.717]; % x, y, z
% View Up Vector
V = [0.692896, -0.378331, -0.613808]; % x, y, z


% first step: rotate and center, such that vertices are aligned with the
% view reference frame

% View Plane Normal 
N = (fp - VRP) / norm(fp-VRP, 2);
% Third Axis
U = -cross(N, V);

% Transformation Matrix
R = [U', V', N'];
t = [-U*VRP', -V*VRP', -N*VRP'];
T = eye(4);
T(1:3, 1:3) = R;
T(4, 1:3) = t;

% Transform points into viewing coordinates
pts_uvn = quickTF(pts_xyz, T);

% pcshow(pts_uvn);
% xlabel('U (Side)');
% ylabel('ViewUp');
% zlabel('Normal');

%% define the view volume
w_l = -50; % left
w_r = 50; % right
w_t = 50; % top
w_b = -50; % bottom

% redefine clipping plane
near = min(pts_uvn(:, 3));
far = max(pts_uvn(:, 3));

% create the normalization matrix
N_mat = eye(4);

N_mat(1, 1) = 2 / (w_r - w_l);
N_mat(2, 2) = 2 / (w_t - w_b);
N_mat(3, 3) = 2 / (far - near); % far - near;

N_mat(4, 1) = - (w_r + w_l) / (w_r - w_l);
N_mat(4, 2) = - (w_t + w_b) / (w_t - w_b);
N_mat(4, 3) = - (far + near) / (far - near);

% normalize all points
pts_norm = quickTF(pts_uvn, N_mat);

% show normalized volume
% pcshow(pts_norm);
% xlabel('U (Side)');
% ylabel('ViewUp');
% zlabel('Normal');

%% Viewport Transformation
% after this, all points are in "soft" pixel coordinates, and have a depth
% value between 0 to 1 assigned to them. 
Nx = 1000; % pixels, starting from 0
Ny = 1000; % pixels

V_mat = eye(4);
V_mat(1, 1) = Nx/2;
V_mat(2, 2) = Ny/2;
V_mat(3, 3) = 1/2;
V_mat(4, 1) = Nx/2;
V_mat(4, 2) = Nx/2;
V_mat(4, 3) = 1/2;

pts_vp = quickTF(pts_norm, V_mat);

% show normalized volume
% pcshow(pts_vp);
% xlabel('U (Side)');
% ylabel('ViewUp');
% zlabel('Normal');

%% Rasterization
% we now iterate over all triangles
ID_map = nan(Nx, Ny); % our "pixels", where each pixel holds an object ID
z_buffer = 2*ones(Nx, Ny);

num_faces = size(fcs, 1);

for t_idx = 1:num_faces % t_idx is ID of triangle
    v_idx = fcs(t_idx, :); % v_idx contains 3 vertex IDs
    p3 = pts_vp(v_idx, :); % p3 is a 3x3 matrix containing one triangle
    p2 = p3(:, 1:2);
    
    % find bounding box in x, y
    xmin = int32(floor(min(p3(:, 1))));
    xmax = int32(ceil(max(p3(:, 1))));
    ymin = int32(floor(min(p3(:, 2))));
    ymax = int32(ceil(max(p3(:, 2))));
    
    % loop over the pixels in the bounding box
    for x = xmin:xmax
        for y = ymin:ymax
            % calculate barycentric coordinates
            a = bary_fn(x, y, p2, 2, 3) / bary_fn(p2(1, 1), p2(1, 2), p2, 2, 3);
            b = bary_fn(x, y, p2, 3, 1) / bary_fn(p2(2, 1), p2(2, 2), p2, 3, 1);
            c = bary_fn(x, y, p2, 1, 2) / bary_fn(p2(3, 1), p2(3, 2), p2, 1, 2);
            
            % if point lies in triangle and in visible field
            if a > 0 && b > 0 && c > 0 && x > 0 && y > 0 && x <= Nx && y <= Ny
                % get depth of point
                z = [a, b, c] * p3(:, 3);
                
                % if depth is smaller than saved depth
                if z < z_buffer(x, y)
                    % buffer new z value
                    z_buffer(x, y) = z;
                    % write object ID to ID_map
                    ID_map(x, y) = t_idx;
                end
            end
        end
    end
    if ~mod(t_idx, 10000)
        fprintf('Finished object %d (%0.2f %%)\n', t_idx, t_idx / num_faces * 100);
    end
end

%% now get a list of all object IDs 
valid_IDs = int32(unique(reshape(ID_map, [], 1))); % nan becomes 0
valid_IDs = valid_IDs(valid_IDs ~= 0);
valid_fcs = fcs(valid_IDs, :);

valid_pt_idx = unique(reshape(valid_fcs, [], 1));
valid_pts = pts_xyz(valid_pt_idx, :);

%% save mesh (Upsampled again) as point cloud

% each pt vector (n x 1) specifies one corner of a triangle for all triangles
pt1 = valid_fcs(:, 1);
pt2 = valid_fcs(:, 2);
pt3 = valid_fcs(:, 3);

% each pos vector (n x 3) holds the coordinates of the pt vector points
pos1 = pts_xyz(pt1, :);
pos2 = pts_xyz(pt2, :);
pos3 = pts_xyz(pt3, :);

% to create a new point on the triangle surface, average the three pos
% vectors
new_pts =  (pos1 + pos2 + pos3) / 3;

%% save the old and new points as new pointcloud
all_pts = single([valid_pts; new_pts]);
clear new_pts pos1 pos2 pos3 pt1 pt2 pt3

pcwrite(pointCloud(all_pts), 'Data/PointClouds/Render.pcd');

%% define barycentric function
function res = bary_fn(x, y, p2, i, j)
    res = (p2(i, 2) - p2(j, 2))*double(x) + ...
          (p2(j, 1) - p2(i, 1)) * double(y) + ...
          p2(i, 1) * p2(j, 2) - p2(j, 1)*p2(i, 2);
end
