function matches = getMatches(descSurface, descModel, par)
% this function provides the matching functionality present in
% GetSphericalDescriptors.m except for Mahalanobis Distance

if isfield(par, 'VERBOSE')
    VERBOSE = par.VERBOSE;
else
    VERBOSE = 1;
end

%% descriptor weighting, normalization, and matching options
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply the normalization and match %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- Un-normalization by adding constant element to end --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.UNNORMALIZE
    % get average descriptor length
    avg_desc_len = mean(vecnorm([descSurface; descModel], 1, 2));
    descSurfaceN = [descSurface,  par.norm_factor*avg_desc_len*ones(size(descSurface, 1), 1)]; 
    descModelN = [descModel,  par.norm_factor*avg_desc_len*ones(size(descModel, 1), 1)];
else
    descSurfaceN = descSurface;
    descModelN = descModel;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- raise descriptors to power before L1-matching --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.CHANGE_METRIC
    descSurfaceC = descSurfaceN.^par.metric_factor;
    descModelC = descModelN.^par.metric_factor;    
else
    descSurfaceC = descSurfaceN;
    descModelC = descModelN;
end

% clear everything that is not needed anymore (intermediate variables)
clear descSurface descModel
clear descSurfaceN descModelN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Match features between Surface and Model / Random Crop %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matches = matchFeatures(descSurfaceC, descModelC, ...
        'Method', par.Method, ...
        'MatchThreshold', par.MatchThreshold, ... 
        'MaxRatio', par.MaxRatio, ... 
        'Metric', par.Metric, ...
        'Unique', par.Unique); 
if VERBOSE
    fprintf('Calculated matches in %0.1f seconds...\n', toc);
end
