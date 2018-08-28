function [T, varargout] = ransac( pts1,pts2,ransacCoef,funcFindTransf,funcDist, REFINE)
    %[T, inlierIdx] = ransac1( pts1,pts2,ransacCoef,funcFindTransf,funcDist )
    %	Use RANdom SAmple Consensus to find a fit from PTS1 to PTS2.
    %	PTS1 is M*n matrix including n points with dim M, PTS2 is N*n;
    %	The transform, T, and the indices of inliers, are returned.
    %
    %	RANSACCOEF is a struct with following fields:
    %	minPtNum,iterNum,thDist,thInlrRatio
    %	MINPTNUM is the minimum number of points with whom can we 
    %	find a fit. For line fitting, it's 2. For homography, it's 4.
    %	ITERNUM is the number of iteration, THDIST is the inlier 
    %	distance threshold and ROUND(THINLRRATIO*n) is the inlier number threshold.
    %
    %	FUNCFINDF is a func handle, f1 = funcFindF(x1,y1)
    %	x1 is M*n1 and y1 is N*n1, n1 >= ransacCoef.minPtNum
    %	f1 can be of any type.
    %	FUNCDIST is a func handle, d = funcDist(f,x1,y1)
    %	It uses f returned by FUNCFINDF, and return the distance
    %	between f and the points, d is 1*n1.
    
    nout = max(nargout, 1);

    minPtNum = ransacCoef.minPtNum;
    iterNum = ransacCoef.iterNum;
    thInlrRatio = ransacCoef.thInlrRatio;
    thDist = ransacCoef.thDist;
    ptNum = size(pts1,1);
    thInlr = round(thInlrRatio*ptNum);

    inlrNum = zeros(1,iterNum);
    inlrNum_refined = zeros(1,iterNum);
    TForms = cell(1,iterNum);

    for p = 1:iterNum 
        % 1. fit using  random points
        randomInts = randperm(ptNum);
        sampleIdx = randomInts(1:minPtNum);
        f1 = funcFindTransf(pts1(sampleIdx, :),pts2(sampleIdx, :));

        % 2. count the inliers, if more than thInlr, refit; else iterate
        dist = funcDist(f1,pts1,pts2);
        inlier1 = find(dist < thDist);
        inlrNum(p) = length(inlier1);

        % refit if enough inliers, check that threshold is still met
        if length(inlier1) >= thInlr
            if REFINE
                f1_ref = funcFindTransf(pts1(inlier1, :),pts2(inlier1, :));
                dist = funcDist(f1_ref,pts1,pts2);
                inlier_refined = find(dist < thDist);
                inlrNum_refined(p) = length(inlier_refined);
                if inlrNum_refined(p) >= thInlr
                    TForms{p} = f1_ref;
                end
            else
                TForms{p} = f1;
            end
        end      
    end

    % 3. choose the coef with the most inliers
    if REFINE
        [maxInliers,idx] = max(inlrNum_refined);
    else
        [maxInliers,idx] = max(inlrNum);
    end
    
    T = TForms{idx};
    FAILED = false;
    try
        dist = funcDist(T,pts1,pts2);
    catch
        % error('RANSAC could not find an appropriate transformation');
        FAILED = true;
        T = [];
        inlierIdx = [];
        numSuccess = 0;
        maxInliers = 0;
    end
    
    if ~FAILED
        inlierIdx = find(dist < thDist);

        if REFINE
            numSuccess = sum(inlrNum_refined >= thInlr);
        else
            numSuccess = sum(inlrNum >= thInlr);
        end

        fprintf('RANSAC succeeded %d times with a maximum of %d Inliers (%0.2f %%)\n', numSuccess, maxInliers, 100*maxInliers/ptNum);
    end
	
    % outputs        
    if nout > 1
        varargout{1} = inlierIdx;
    end
    if nout > 2
        varargout{2} = numSuccess;
    end 
    if nout > 3
        varargout{3} = maxInliers;
    end 
    if nout > 4
        varargout{4} = 100*maxInliers/ptNum;
    end 
    
end