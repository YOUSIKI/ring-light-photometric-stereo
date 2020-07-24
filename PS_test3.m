%3ͬ�ȼ�����շ��򣬷ֲ����Գƣ����跽��ƽ�У�ֻ�о���˥��-��ƽ��
clear; close all; clc;
addpath(genpath(pwd));  % ��ӵ�ǰ·���µ�������Ŀ¼

g_pic_num = 3;
shadowThresh = 0.01;

%% ���շ���
L = zeros(3, g_pic_num);
Slant = 45;
Slant_sin = sin(Slant/180*pi);
Slant_cos = cos(Slant/180*pi);
for i = 1:g_pic_num
    Tilt = (i-1) * 360/g_pic_num;
    %%��ʵ�������ϵ
    L(1,i) = Slant_sin * cos(Tilt/180*pi);
    L(2,i) = Slant_sin * sin(Tilt/180*pi);
    L(3,i) = Slant_cos;
end

%% ͼ�����,����һ��ƽ�棬����0��0��1��
w1 = 600; h1 = 400;
[X,Y]= meshgrid(1:w1, 1:h1);

x = (X - 300)/20;
y = (Y - 200)/20;

r = 40; h = 40;
d1 = sqrt( (x - r).^2 + y.^2 + h^2 );
d2 = sqrt( (x + r/2).^2 + (y-sqrt(3)*r/2).^2 + h^2);
d3 = sqrt( (x + r/2).^2 + (y+sqrt(3)*r/2).^2 + h^2);

I1 = 255./(d1.^2);
I2 = 255./(d2.^2);
I3 = 255./(d3.^2);
figure; subplot(1,3,1); imshow(I1); subplot(1,3,2); imshow(I2); subplot(1,3,3); imshow(I3);

%% ����ͼƬ
[M, N, C] = size(I1);
I = ones(M, N, g_pic_num);
% max_I = max([max(max(max(I1))),max(max(max(I2))),max(max(max(I3)))]); 
max_I = max([prctile(I1(:), 99),prctile(I2(:), 99),prctile(I3(:), 99)]);
I(:,:,1) = I1/max_I;
I(:,:,2) = I2/max_I;
I(:,:,3) = I3/max_I;

%% ͼ��ĳ߶���Ϣ
[g_rows, g_cols] = size(I1);
g_length = g_rows * g_cols;

% Create a shadow mask.
shadow_mask = (I > shadowThresh);
se = strel('disk', 2);
for i = 1:g_pic_num
  % Erode the shadow map to handle de-mosaiking artifact near shadow boundary.
  shadow_mask(:,:,i) = imerode(shadow_mask(:,:,i), se);
end

[rho, n] = PhotometricStereo(I, shadow_mask, L);

%% Visualize the normal map. axis xy;
N_RGB(:,:,1)  = (n(:,:,1) + 1) / 2;
N_RGB(:,:,2)  = (n(:,:,2) + 1) / 2;
N_RGB(:,:,3)  = n(:,:,3);
figure; imshow(N_RGB); 
% figure; imshow(rho);

%% Estimate depth map from the normal vectors.
fprintf('Estimating depth map from normal vectors...\n');
p = -n(:,:,1) ./ n(:,:,3);
q = -n(:,:,2) ./ n(:,:,3);
p(isnan(p)) = 0; %�ж������Ԫ���Ƿ���NaN�� NaN �� Not a Number ����д��
q(isnan(q)) = 0;

figure; subplot(1,2,1); mesh(p); subplot(1,2,2); mesh(q);

%% integration
Height_poisson =poisson_solver_function_neumann(p, q);
% Height_poisson = flipud(Height_poisson);
figure;mesh(Height_poisson);

%% ��������
% ��� y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;
index = floor(randi(length(Height_poisson(:)),100,1));
Height_diff = Height_poisson(index) - 0;
[y, x] = ind2sub(size(Height_poisson), index); %ind2sub�õ��������У�ע�� ���� �� xy �෴
AA = [x.^2, y.^2, x, y, ones(length(x),1)];

% ����1���������
para = AA \ Height_diff;

% ���ݲ���������ʵ���
[m, n] = size(Height_poisson); %size �õ���������
[xx, yy] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ
fitted_height = [xx(:).^2, yy(:).^2, xx(:), yy(:), ones(m*n,1)] * para(:);
Height_correct = Height_poisson(:) - fitted_height;
Height_correct = reshape(Height_correct, m, n);
fitted_height = reshape(fitted_height, m, n);
figure; mesh(fitted_height);
figure; mesh(Height_correct);

% ������
Height_d = Height_correct - 0;
figure; mesh(Height_d); title('������');