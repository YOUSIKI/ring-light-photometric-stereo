%   Title: A new ring-light photometric stereo compensation method
%
%   Author: Hao Fan.
%   Created: March 25, 2020.

clear; close all; clc;
addpath(genpath(pwd));

%% Simulation for a flat plane. normal=(0,0,1) (图像仿真,仿真一个平面)
w1 = 600; h1 = 600;
[X,Y]= meshgrid(1:w1, 1:h1);

scale = 3; %6-5,3-10,2-15,1.5-20
x = (X - 300)/scale; %+-15 length(长)
y = (Y - 300)/scale; %+-15 width(宽)

r = 400; h = 400;
d0 = sqrt( (x - r).^2 + y.^2 + h^2 );
d30 = sqrt( (x - sqrt(3)*r/2).^2 + (y - r/2).^2 + h^2 );
d60 = sqrt( (x - r/2).^2 + (y-sqrt(3)*r/2).^2 + h^2 );
d90 = sqrt( (x).^2 + (y - r).^2 + h^2 );
d120 = sqrt( (x + r/2).^2 + (y-sqrt(3)*r/2).^2 + h^2 );
d150 = sqrt( (x + sqrt(3)*r/2).^2 + (y - r/2).^2 + h^2 );
d180 = sqrt( (x + r).^2 + y.^2 + h^2 );
d210 = sqrt( (x + sqrt(3)*r/2).^2 + (y+r/2).^2 + h^2 );
d240 = sqrt( (x + r/2).^2 + (y+sqrt(3)*r/2).^2 + h^2 );
d270 = sqrt( (x).^2 + (y+r).^2 + h^2 );
d300 = sqrt( (x - r/2).^2 + (y+sqrt(3)*r/2).^2 + h^2 );
d330 = sqrt( (x - sqrt(3)*r/2).^2 + (y+r/2).^2 + h^2 );

I_base = 10000; para_attenuation = 3;
I0 = I_base./(d0.^para_attenuation);
I60 = I_base./(d60.^para_attenuation);
I120 = I_base./(d120.^para_attenuation);
I180 = I_base./(d180.^para_attenuation);
I240 = I_base./(d240.^para_attenuation);
I300 = I_base./(d300.^para_attenuation);
I30 = I_base./(d30.^para_attenuation);
I90 = I_base./(d90.^para_attenuation);
I150 = I_base./(d150.^para_attenuation);
I210 = I_base./(d210.^para_attenuation);
I270 = I_base./(d270.^para_attenuation);
I330 = I_base./(d330.^para_attenuation);

% Adjust the imaging value in [0 255] (图像亮度调整到[0 255])
max_I1 = max([prctile(I0(:), 99),prctile(I60(:), 99),prctile(I120(:), 99),prctile(I180(:), 99),prctile(I240(:), 99),prctile(I300(:), 99)]);
max_I2 = max([prctile(I30(:), 99),prctile(I90(:), 99),prctile(I150(:), 99),prctile(I210(:), 99),prctile(I270(:), 99),prctile(I330(:), 99)]);
max_I = max(max_I1, max_I2);

I0 = I0/max_I * 255;
I60 = I60/max_I * 255;
I120 = I120/max_I * 255;
I180 = I180/max_I * 255;
I240 = I240/max_I * 255;
I300 = I300/max_I * 255;
I30 = I30/max_I * 255;
I90 = I90/max_I * 255;
I150 = I150/max_I * 255;
I210 = I210/max_I * 255;
I270 = I270/max_I * 255;
I330 = I330/max_I * 255;
figure; subplot(2,3,1); imshow(uint8(I0)); title('0°'); subplot(2,3,2); imshow(uint8(I60)); title('60°'); subplot(2,3,3); imshow(uint8(I120)); title('120°');
subplot(2,3,4); imshow(uint8(I180)); title('180°'); subplot(2,3,5); imshow(uint8(I240)); title('240°'); subplot(2,3,6); imshow(uint8(I300)); title('300°');

%% ring-light photometric stereo with 6 lights following Quadratic attenuation
% 6对称光照方向，逆平方距离衰减
shadowThresh = 0.01;
[M, N, C] = size(I0);

m_symmetry = 0; %1对称，0不对称asymmetry
g_pic_num = 9;

%% 12对称图，6对称图，4对称图, 3非对
if m_symmetry == 1
    if g_pic_num == 12
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I30;    I(:,:,3) = I60;
        I(:,:,4) = I90;    I(:,:,5) = I120;    I(:,:,6) = I150;
        I(:,:,7) = I180;    I(:,:,8) = I210;    I(:,:,9) = I240;
        I(:,:,10) = I270;    I(:,:,11) = I300;    I(:,:,12) = I330;
    elseif g_pic_num == 6  %6个对称
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I60;    I(:,:,3) = I120;
        I(:,:,4) = I180;    I(:,:,5) = I240;    I(:,:,6) = I300;
    elseif g_pic_num == 4  %4个对称
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I90;
        I(:,:,3) = I180;    I(:,:,4) = I270;
    elseif g_pic_num == 3  %3个对称
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I120;
        I(:,:,3) = I240; 
    end
    
    %% main light direction (光照方向)
    L = zeros(3, g_pic_num);
    Slant = 45;
    Slant_sin = sin(Slant/180*pi);
    Slant_cos = cos(Slant/180*pi);
    for i = 1:g_pic_num
        Tilt = (i-1) * 360/g_pic_num;
        
        L(1,i) = Slant_sin * cos(Tilt/180*pi);
        L(2,i) = Slant_sin * sin(Tilt/180*pi);
        L(3,i) = Slant_cos;
    end
    
elseif m_symmetry == 0
    % main light direction (光照方向)
    L_base = zeros(3, 12);
    Slant = 45;
    Slant_sin = sin(Slant/180*pi);
    Slant_cos = cos(Slant/180*pi);
    for i = 1:12
        Tilt = (i-1) * 360/12;
        
        L_base(1,i) = Slant_sin * cos(Tilt/180*pi);
        L_base(2,i) = Slant_sin * sin(Tilt/180*pi);
        L_base(3,i) = Slant_cos;
    end
    %%%——————————————————————————————%%%
    if g_pic_num == 9 %1-9/12
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I30;    I(:,:,3) = I60;
        I(:,:,4) = I90;    I(:,:,5) = I120;    I(:,:,6) = I150;
        I(:,:,7) = I180;    I(:,:,8) = I210;    I(:,:,9) = I240;
        L = L_base(:,1:9);
    end
    if g_pic_num == 6 %6不对称
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I90;    I(:,:,3) = I120;
        I(:,:,4) = I180;    I(:,:,5) = I240;    I(:,:,6) = I330;
        
        L(:,1) = L_base(:,1);  L(:,2) = L_base(:,4);  L(:,3) = L_base(:,5);
        L(:,4) = L_base(:,7);  L(:,5) = L_base(:,9);  L(:,6) = L_base(:,12);
    end
    if g_pic_num == 4 %4不对称
        I = ones(M, N, g_pic_num);
        I(:,:,1) = I0;    I(:,:,2) = I30;    I(:,:,3) = I180;
        I(:,:,4) = I300;   
        
        L(:,1) = L_base(:,1);  L(:,2) = L_base(:,2);  L(:,3) = L_base(:,7);
        L(:,4) = L_base(:,11);
    end
else
    return;
end

%% Size of the image (图像的尺度信息)
[g_rows, g_cols] = size(I0);
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
figure; imshow(N_RGB); title('normal');
% figure; imshow(rho);

%% Estimate depth map from the normal vectors.
fprintf('Estimating depth map from normal vectors...\n');
p = -n(:,:,1) ./ n(:,:,3);
q = -n(:,:,2) ./ n(:,:,3);
p(isnan(p)) = 0;
q(isnan(q)) = 0;

figure; subplot(1,2,1); mesh(p); title('Gradient p'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('p(pixel)');
subplot(1,2,2); mesh(q); title('Gradient q'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('q(pixel)');
set(gcf,'unit','centimeters','position',[10 5 13 6]);

%% integration
Height_poisson =poisson_solver_function_neumann(p, q);
% Height_poisson = flipud(Height_poisson);
% figure;mesh(Height_poisson); title('积分高度(像素)')

%% The propsed method to fit the bias. (拟合 y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;)
index = floor(randi(length(Height_poisson(:)),500,1));
Height_diff = Height_poisson(index) - 0; % 像素尺度
[y, x] = ind2sub(size(Height_poisson), index); %ind2sub得到的是行列，注意 行列 与 xy 相反
AA = [x.^2, y.^2, x, y, ones(length(x),1)];

% 方法1：矩阵除法
para = AA \ Height_diff;

% Compute the corrected height (根据参数计算真实结果)
[m, n] = size(Height_poisson); %size 得到的是行列
[xx, yy] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系 % 注意像素尺度，不能混淆！！
fitted_height = [xx(:).^2, yy(:).^2, xx(:), yy(:), ones(m*n,1)] * para(:);
Height_correct = Height_poisson(:) - fitted_height;
Height_correct = reshape(Height_correct, m, n);
fitted_height = reshape(fitted_height, m, n);

figure; mesh(Height_correct); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Corrected Height(pixel)'); grid on; % axis equal;
figure; mesh(fitted_height); title('Fitted Bias Height(pixel)');

figure; subplot(1,2,1); mesh(Height_poisson); title('Initial Height'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Initial Height(pixel)');
subplot(1,2,2); mesh(Height_correct); title('Correctted Height'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Corrected Height(pixel)');
set(gcf,'unit','centimeters','position',[10 5 13 6]);

% Height Error Analysis
Height_d = Height_correct - 0; % 像素尺度
disp(['平均高度误差：' num2str(mean(mean(abs(Height_d)))) '像素']);

%figure; mesh(Height_d/scale); title('Height error(mm)');colorbar;
disp(['平均高度误差：' num2str(mean(mean(abs(Height_d/scale)))) '毫米']);

% Angle error Analysis (偏差角度的误差分析)
% 法向没有尺度，但是要确认光度立体高度积分的尺度都是像素，或都是真实尺度；混淆后计算的法向是错误的！！！
[p_correct,q_correct] = gradient(Height_correct);
bb = sqrt(p_correct.^2 + q_correct.^2 + 1);
% Ninit(:,:,1) = - p_correct./bb;
% Ninit(:,:,2) = - q_correct./bb;
% Ninit(:,:,3) = 1./bb;
Height_n_error = acos(1./bb)/pi * 180;

figure; subplot(1,2,1); mesh(Height_d); title('Height error'); xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Height error(pixel)');
subplot(1,2,2);mesh(Height_n_error); title('Angle error');  xlabel('x(pixel)'); ylabel('y(pixel)'); zlabel('Angle error(pixel)');
set(gcf,'unit','centimeters','position',[10 5 13 6]);

disp(['平均角度误差：' num2str(mean(mean(Height_n_error))) '度']);