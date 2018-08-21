function pcOut = pcRigidBodyTF(pcIn, r, t)
% performs a 3D rigid body transform given the 
% rotation r = [rx, ry, rz] in rad and the
% translation t = [tx, ty, tz] in the units used by pcIn

    %% define 4x4 affine transformation
    if isequal(size(r), [3, 1])
        r = r';
    elseif ~isequal(size(r), [1, 3])
        error('Rotation must be either 1x3 or 3x1 matrix.');
    end
        
    R = eul2rotm(r, 'XYZ');
    T = eye(4);
    T(1:3, 1:3) = R;
    T(4, 1:3) = t;
    tform = affine3d(T);

    pcOut = pctransform(pcIn, tform);
end