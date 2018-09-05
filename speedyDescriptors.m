%% function to calculate a descriptor for each point
function [feat, desc] = speedyDescriptors(pts, spts, options)
    % speedy implementation around getSpacialHistogramDescriptorsm.
    % using division into refions to speed up the calculation process
    
    % returns:
        % - feat: feature locations
        % - desc: feature descriptors        

    max_region_size = options.max_region_size;    
    VERBOSE = options.VERBOSE;
    options.VERBOSE = options.VERBOSE - 1;
    
    % divide Model into regions
    pc = pointCloud(pts);
    xyzLim = [pc.XLimits; pc.YLimits; pc.ZLimits];
    xyzRange = [xyzLim(1, 2) - xyzLim(1, 1); xyzLim(2, 2) - xyzLim(2, 1); xyzLim(3, 2) - xyzLim(3, 1)];
    xyzNumRegions = ceil(xyzRange / max_region_size);
    xyzRegionSize = xyzRange ./ xyzNumRegions;

    % boundary grid in all three axes for the cuboid regions
    xBounds = xyzLim(1, 1):xyzRegionSize(1):xyzLim(1, 2);
    yBounds = xyzLim(2, 1):xyzRegionSize(2):xyzLim(2, 2);
    zBounds = xyzLim(3, 1):xyzRegionSize(3):xyzLim(3, 2);
    assert(isequal([size(xBounds, 2),size(yBounds, 2),size(zBounds, 2)]'-1, xyzNumRegions));


    feat = [];
    desc = [];
    R = options.R;
    
    total_num_regions = prod(xyzNumRegions);
    if VERBOSE
        fprintf('Divided model into %d regions\n', total_num_regions);
    end
    
    % create waitbar
    waitbar_handle = waitbar(0,'');
    progress = 0;
    
    tic
    for ix = 1:xyzNumRegions(1)
        for iy = 1:xyzNumRegions(2)
            for iz = 1:xyzNumRegions(3)
                % get points in the respective region with margin of descOptM.R = 3.5;
                mask = pts(:, 1) > xBounds(ix)-R & pts(:, 1) < xBounds(ix+1)+R ...
                     & pts(:, 2) > yBounds(iy)-R & pts(:, 2) < yBounds(iy+1)+R...
                     & pts(:, 3) > zBounds(iz)-R & pts(:, 3) < zBounds(iz+1)+R;
                mask = cat(2, mask, mask, mask);
                pts_crop = reshape(pts(mask), [], 3);    

                % same for the sample points, but without the margin of R
                mask = spts(:, 1) > xBounds(ix)-R & spts(:, 1) < xBounds(ix+1)+R ...
                 & spts(:, 2) > yBounds(iy)-R & spts(:, 2) < yBounds(iy+1)+R...
                 & spts(:, 3) > zBounds(iz)-R & spts(:, 3) < zBounds(iz+1)+R;
                mask = cat(2, mask, mask, mask);
                spts_crop = reshape(spts(mask), [], 3); 

                % use only those points for descriptor calculation
                [feat_crop, desc_crop] = ...
                        getSpacialHistogramDescriptors(pts_crop, spts_crop, options);

                % put the results together
                feat = [feat; feat_crop];
                desc = [desc; desc_crop];
                
                % update waitbar
                progress = progress + 1;
                progress_frac = progress/total_num_regions;
                t_total = double(toc)/progress_frac;
                t_remain = round(t_total - toc);
                msg = sprintf('Progress: %0.2f %%, time remaining: %d seconds', progress_frac*100, t_remain); 
                waitbar(progress_frac,waitbar_handle,msg)
            end
        end
    end
    if VERBOSE
        fprintf('Calculated descriptors in %0.1f seconds...\n', toc);
    end
    close(waitbar_handle);
end
