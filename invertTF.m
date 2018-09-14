function TFinv = invertTF(TF)
% this function inverts a  4x4 transform for 3D points such that 
% TFinv = inv(TF);

TFinv = eye(4);
TFinv(1:3, 1:3) = TF(1:3, 1:3)';
TFinv(4, 1:3) = -TF(4, 1:3) * TF(1:3, 1:3)';
end