%6�Գƹ��շ��򣬼��跽��ƽ�У�ֻ�о���˥��-��ƽ��
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
w1 = 600; h1 = 600;
[X,Y]= meshgrid(1:w1, 1:h1);

scale = 3; %6-5,3-10,2-15,1.5-20
x = (X - 300)/scale; %+-15��
y = (Y - 300)/scale; %+-15��

r = 400; h = 400;
d1 = sqrt( (x - r).^2 + y.^2 + h^2 );
d2 = sqrt( (x - r/2).^2 + (y-sqrt(3)*r/2).^2 + h^2 );
d3 = sqrt( (x + r/2).^2 + (y-sqrt(3)*r/2).^2 + h^2 );
d4 = sqrt( (x + r).^2 + y.^2 + h^2 );
d5 = sqrt( (x + r/2).^2 + (y+sqrt(3)*r/2).^2 + h^2 );
d6 = sqrt( (x - r/2).^2 + (y+sqrt(3)*r/2).^2 + h^2 );

I0 = 10000; para_attenuation = 3;
I1 = I0./(d1.^para_attenuation);
I2 = I0./(d2.^para_attenuation);
I3 = I0./(d3.^para_attenuation);
I4 = I0./(d4.^para_attenuation);
I5 = I0./(d5.^para_attenuation);
I6 = I0./(d6.^para_attenuation);

% scale1 = exp(-0.1*d1) ./ ((d1).^3);
% I0 = 500000;
% I1 = I0.* exp(-0.1*d1)./ ((d1).^3);
% I2 = I0.* exp(-0.1*d2)./ ((d2).^3);
% I3 = I0.* exp(-0.1*d3)./ ((d3).^3);
% I4 = I0.* exp(-0.1*d4)./ ((d4).^3);
% I5 = I0.* exp(-0.1*d5)./ ((d5).^3);
% I6 = I0.* exp(-0.1*d6)./ ((d6).^3);

% ��ʾ��������[0 255]
max_I = max([prctile(I1(:), 99),prctile(I2(:), 99),prctile(I3(:), 99),prctile(I4(:), 99),prctile(I5(:), 99),prctile(I6(:), 99)]);
I1 = I1/max_I * 255;
I2 = I2/max_I * 255;
I3 = I3/max_I * 255;
I4 = I4/max_I * 255;
I5 = I5/max_I * 255;
I6 = I6/max_I * 255;
figure; subplot(2,3,1); imshow(uint8(I1)); title('0��'); subplot(2,3,2); imshow(uint8(I2)); title('60��'); subplot(2,3,3); imshow(uint8(I3)); title('120��');
subplot(2,3,4); imshow(uint8(I4)); title('180��'); subplot(2,3,5); imshow(uint8(I5)); title('240��'); subplot(2,3,6); imshow(uint8(I6)); title('300��');

%% ����ͼƬ
[M, N, C] = size(I1);
I = ones(M, N, g_pic_num);
I(:,:,1) = I1;
I(:,:,2) = I2;
I(:,:,3) = I3;
I(:,:,4) = I4;
I(:,:,5) = I5;
I(:,:,6) = I6;
%6���Գƣ�5�����Գ�
g_pic_num_2 = 6;
I_new = I;
L_new = L;

% % 4���Գ�, 3�����Գ�
% g_pic_num_2 = 3;
% I_new = ones(M, N, g_pic_num_2);
% I_new(:,:,1) = I1;
% I_new(:,:,2) = I3;
% I_new(:,:,3) = I5;
% % I_new(:,:,4) = I4;
% L_new = zeros(3, g_pic_num_2);
% L_new(:,1) = L(:,1);
% L_new(:,2) = L(:,3);
% L_new(:,3) = L(:,5);
% % L_new(:,4) = L(:,4);

%% ͼ��ĳ߶���Ϣ
[g_rows, g_cols] = size(I1);
g_length = g_rows * g_cols;

% Create a shadow mask.
shadow_mask = (I_new > shadowThresh);
se = strel('disk', 2);
for i = 1:g_pic_num_2
  % Erode the shadow map to handle de-mosaiking artifact near shadow boundary.
  shadow_mask(:,:,i) = imerode(shadow_mask(:,:,i), se);
end

[rho, n] = PhotometricStereo(I_new, shadow_mask, L_new);

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

figure; subplot(1,2,1); mesh(p); title('p'); subplot(1,2,2); mesh(q); title('q');

%% integration
Height_poisson =poisson_solver_function_neumann(p, q);
% Height_poisson = flipud(Height_poisson);
% figure;mesh(Height_poisson); title('���ָ߶�(����)')

%% ��������
% ��� y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;
index = floor(randi(length(Height_poisson(:)),500,1));
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
figure; mesh(fitted_height); title('���߶�(����)')

figure; subplot(1,2,1); mesh(Height_poisson); title('Initial Height(pixel)'); subplot(1,2,2); mesh(Height_correct); title('Correct Height(pixel)');

% ������
Height_d = Height_correct - 0;
% figure; mesh(Height_d); title('������(����)');colorbar;% axis equal;
disp(['ƽ���߶���' num2str(mean(mean(abs(Height_d)))) '����']);

figure; mesh(Height_d/scale); title('������(����)');colorbar;
disp(['ƽ���߶���' num2str(mean(mean(abs(Height_d/scale)))) '����']);

% ƫ��Ƕȵ�������
[p_correct,q_correct] = gradient(Height_correct);
bb = sqrt(p_correct.^2 + q_correct.^2 + 1);
% Ninit(:,:,1) = - p_correct./bb;
% Ninit(:,:,2) = - q_correct./bb;
% Ninit(:,:,3) = 1./bb;
Height_n_error = acos(1./bb)/pi * 180;

figure; subplot(1,2,1); mesh(Height_d); title('Height error (pixel)');colorbar;
subplot(1,2,2);mesh(Height_n_error); title('Angle error (degree)'); colorbar;

disp(['ƽ���Ƕ���' num2str(mean(mean(Height_n_error))) '��']);