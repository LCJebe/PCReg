function [T, inlierIdx] = ransac( pts1,pts2,ransacCoef,funcFindTransf,funcDist )
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
    %	For line fitting, it should calculate the dist between the line and the
    %	points [x1;y1]; for homography, it should project x1 to y2 then
    %	calculate the dist between y1 and y2.
    %	Yan Ke @ THUEE, 20110123, xjed09@gmail.com


    minPtNum = ransacCoef.minPtNum;
    iterNum = ransacCoef.iterNum;
    thInlrRatio = ransacCoef.thInlrRatio;
    thDist = ransacCoef.thDist;
    ptNum = size(pts1,1);
    thInlr = round(thInlrRatio*ptNum);

    inlrNum = zeros(1,iterNum);
    TForms = cell(1,iterNum);

    for p = 1:iterNum %PARFOR
        % 1. fit using  random points
        randomInts = randperm(ptNum);
        sampleIdx = randomInts(1:minPtNum);
        f1 = funcFindTransf(pts1(sampleIdx, :),pts2(sampleIdx, :));

        % 2. count the inliers, if more than thInlr, refit; else iterate
        dist = funcDist(f1,pts1,pts2);
        inlier1 = find(dist < thDist);
        inlrNum(p) = length(inlier1);

        % refit if enough inliers
        if length(inlier1) >= thInlr
            TForms{p} = funcFindTransf(pts1(inlier1, :),pts2(inlier1, :));
        end      
    end

    % 3. choose the coef with the most inliers
    [maxInliers,idx] = max(inlrNum);
    T = TForms{idx};
    try
        dist = funcDist(T,pts1,pts2);
        inlier_refined = find(dist < thDist);
        maxInliers_refined = length(inlier_refined);
    catch
        error('RANSAC could not find an appropriate transformation');
    end
    inlierIdx = find(dist < thDist);
    
    numSuccess = sum(inlrNum >= thInlr);
    
    fprintf('RANSAC succeeded %d times with a maximum of %d Inliers\n', numSuccess, maxInliers_refined);
	
end