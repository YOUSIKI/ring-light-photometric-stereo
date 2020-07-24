%6�Գƹ��շ��򣬷�����ƽ�У���ƽ��-б��
clear; close all; clc;
addpath(genpath(pwd));  % ��ӵ�ǰ·���µ�������Ŀ¼

g_pic_num = 6;
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

x = (X - 300)/20; %+-15��
y = (Y - 200)/20; %+-10��

% б��
z_a = -0.1; z_b = -0.1;
z = -z_a - z_b*y;
figure;mesh(z);title('z');

r = 40; h = 40;
% ��Դ������ľ���仯
d1 = sqrt( (x - r).^2 + y.^2 + (h-z).^2 );
d2 = sqrt( (x - r/2).^2 + (y-sqrt(3)*r/2).^2 + (h-z).^2 );
d3 = sqrt( (x + r/2).^2 + (y-sqrt(3)*r/2).^2 + (h-z).^2 );
d4 = sqrt( (x + r).^2 + y.^2 + (h-z).^2 );
d5 = sqrt( (x + r/2).^2 + (y+sqrt(3)*r/2).^2 + (h-z).^2 );
d6 = sqrt( (x - r/2).^2 + (y+sqrt(3)*r/2).^2 + (h-z).^2 );

I0 = 10000;
% б��ķ����ǿ N*L/|L|,����|N|��ͬδ����
I1 = I0*(z_a .*(r-x) + z_b.*(-y)+(h-z))./(d1);
I2 = I0*(z_a .*(r/2-x) + z_b.*(sqrt(3)*r/2-y)+(h-z))./(d2);
I3 = I0*(z_a .*(-r/2-x) + z_b.*(sqrt(3)*r/2-y)+(h-z))./(d3);
I4 = I0*(z_a .*(-r-x) + z_b.*(-y)+(h-z))./(d4);
I5 = I0*(z_a .*(-r/2-x) + z_b.*(-sqrt(3)*r/2-y)+(h-z))./(d5);
I6 = I0*(z_a .*(r/2-x) + z_b.*(-sqrt(3)*r/2-y)+(h-z))./(d6);

% ��ʾ��������[0 255]
max_I = max([prctile(I1(:), 99),prctile(I2(:), 99),prctile(I3(:), 99),prctile(I4(:), 99),prctile(I5(:), 99),prctile(I6(:), 99)]);
I1 = I1/max_I * 255;
I2 = I2/max_I * 255;
I3 = I3/max_I * 255;
I4 = I4/max_I * 255;
I5 = I5/max_I * 255;
I6 = I6/max_I * 255;
figure; subplot(2,3,1); imshow(uint8(I1)); subplot(2,3,2); imshow(uint8(I2)); subplot(2,3,3); imshow(uint8(I3));
figure; subplot(2,3,4); imshow(uint8(I4)); subplot(2,3,5); imshow(uint8(I5)); subplot(2,3,6); imshow(uint8(I6));

%% ����ͼƬ
[M, N, C] = size(I1);
I = ones(M, N, g_pic_num);
I(:,:,1) = I1;
I(:,:,2) = I2;
I(:,:,3) = I3;
I(:,:,4) = I4;
I(:,:,5) = I5;
I(:,:,6) = I6;

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
index = floor(randi(length(Height_poisson(:)),1000,1));
Height_diff = Height_poisson(index) - z(index);  %��ͬ��
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
Height_d = Height_correct - z;
figure; mesh(Height_d); title('������');