%% load colored and uncolored model
pc_col = pcread('Data/PointClouds/Model.pcd');
pc = pcread('Data/PointClouds/ModelSmoothUp3.pcd');


%% find nearest neighbor for each point and transfer colors

colors = pc_col.Color;
new_colors = zeros(pc.Count, 3);

% create waitbar
f = waitbar(0,'');

tic
for i = 1:pc.Count
    point = pc.Location(i,:);
    [idx,~] = findNearestNeighbors(pc_col,point,1);
    new_colors(i, :) = colors(idx, :);
    % update waitbar
    if mod(i, 1000) == 0
        progress = double(i)/pc.Count*100;
        t_total = toc*pc.Count/i;
        t_remain = round(t_total - toc);
        msg = sprintf('Progress: %0.2f %%, time remaining: %d seconds', progress, t_remain); 
        waitbar(progress/100,f,msg)
    end
end
close(f);

%% save centered pointcloud with new colors
pc_new = pointCloud(pc.Location, 'Color', uint8(new_colors));

pcwrite(pc_new, 'Data/PointClouds/ModelSmoothColorUp3.pcd');
