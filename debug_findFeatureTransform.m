% load features with random transform, and from raw model

f_rand = load('Data/Descriptors/featSurface0.3_rand.mat');
f_rand = f_rand.featSurface;
num_rand = size(f_rand, 1);
f_raw = load('Data/Descriptors/featSurface0.3_raw.mat');
f_raw = f_raw.featSurface;
num_raw = size(f_raw, 1);

% remove mean (centroid)
f_rand_c = f_rand - mean(f_rand, 1);
f_raw_c = f_raw - mean(f_raw, 1);

% perform pca on both
[coeff_rand, lrf_rand, ~] = pca(f_rand_c, 'Algorithm', 'eig');
[coeff_raw, lrf_raw, ~] = pca(f_raw_c, 'Algorithm', 'eig');

% sign issues: rand
k = num_rand;
x_sign = sum(sign(lrf_rand(:, 1)) == 1) >= k/2;
z_sign = sum(sign(lrf_rand(:, 3)) == 1) >= k/2;
x_sign = x_sign*2-1;
z_sign = z_sign*2-1;
y_sign = det(coeff_rand .* [x_sign, 1, z_sign]);
coeff_rand_u = coeff_rand .* [x_sign, y_sign, z_sign];

% sign issues: raw
k = num_raw;
x_sign = sum(sign(lrf_raw(:, 1)) == 1) >= k/2;
z_sign = sum(sign(lrf_raw(:, 3)) == 1) >= k/2;
x_sign = x_sign*2-1;
z_sign = z_sign*2-1;
y_sign = det(coeff_raw .* [x_sign, 1, z_sign]);
coeff_raw_u = coeff_raw .* [x_sign, y_sign, z_sign];

% transform
f_rand_tf = f_rand_c * coeff_rand_u;
f_raw_tf = f_raw_c * coeff_raw_u;

% display
col1 = ones(num_rand, 3) .* [255, 0, 0];
col2 = ones(num_raw, 3) .* [0, 0, 0];
pcCombined = pointCloud([f_rand_tf; f_raw_tf], 'Color', [col1; col2]);

pcshow(pcCombined);

% get transformation needed from f_rand to f_raw
R = coeff_rand_u * coeff_raw_u';
t = - mean(f_rand, 1) * R +mean(f_raw, 1);

TF = eye(4);
TF(1:3, 1:3) = R;
TF (4, 1:3) = t;

% if everything is right, then [f_rand 1] * TF ==> [f_raw 1]
f_rand_TF = quickTF(f_rand, TF);

% display to prove
pcCombined2 = pointCloud([f_rand_TF; f_raw], 'Color', [col1; col2]);

figure();
pcshow(pcCombined2);

% test the inverse transform
f_raw_TF = quickTF(f_raw, invertTF(TF));

% display to prove
pcCombined3 = pointCloud([f_rand; f_raw_TF], 'Color', [col1; col2]);

figure();
pcshow(pcCombined3);