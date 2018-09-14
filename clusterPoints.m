function clusters = clusterPoints(pts, r)
% function that returns point clusters based on their euclidean distance. 
% all points which distance is < r will end up in one region / cluster

num_pts = size(pts, 1);

% initialize frontier, unexplored, explored (indices into "pts"!)
frontier = 1;
unexplored = 2:num_pts;
explored = [];
num_explored = 0;

num_clusters = 0;
clusters = {};

while num_explored < num_pts
    while ~isempty(frontier)
        % do range search for every point in frontier
        idx = rangesearch(pts, pts(frontier, :), r);
        
        % add frontier to explored
        explored = union(explored, frontier);
        
        % get unique list of found points, add those to for new frontier,
        % if they haven't been explored already
        in_range = unique([idx{:}]);
        frontier = setdiff(in_range, explored);
        
        % remove new frontier from unexplored points
        unexplored = setdiff(unexplored, frontier);
    end
    
    % add all found and explored points to one cluster
    num_clusters = num_clusters + 1;
    clusters{num_clusters} = explored;
    num_explored = num_explored + length(explored);
    
    % the first frontier is empty now. start a new frontier
    % also clear explored
    if ~isempty(unexplored)
        frontier = unexplored(1);
        unexplored = unexplored(2:end);
        explored = [];
    end
end

end



