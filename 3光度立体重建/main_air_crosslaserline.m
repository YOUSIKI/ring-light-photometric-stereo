%%交叉线激光矫正，附带尺度矫正
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   Modify: 2021.4.29
clear; close all; clc;
addpath(genpath(pwd));  % 添加当前路径下的所有子目录

%% 初始化
src = '..\1相机及激光面标定\calibration20210425-paper\';
g_src = '..\0实验数据\20210425-paper\4texture_stone2\'; 
%1cushion;2shell;3texture_stone1;4texture_stone2;5throw_pillow1;6throw_pillow2;7word_stone;8sailboat

load([src 'cameraParams']);
K = (cameraParams.IntrinsicMatrix)';
control_part = 1;  %0为全部，1为局部

image = imread([g_src '0.bmp']); img = image(:,:,1);
ROI = 6; sizeImg = floor(size(img)/ROI);
rect = [sizeImg(2)+1, sizeImg(1)+1, sizeImg(2)*4-1, sizeImg(1)*4-1]; %imcrop(img, rect); 按照xy轴裁剪
rectImg = imcrop(img, rect);
clear img image

% 根据有效区，修改相机参数
K(1,3) = K(1,3) - sizeImg(2)*1; %x -- col
K(2,3) = K(2,3) - sizeImg(1)*1; %y -- row
clear sizeImg

%% -- step 1 -- 获取激光点真实三维数据
disp('step 1 -- 交叉激光线获取');
load ([g_src 'LaserHeightV.mat'], 'LaserHeightV');%laser1 是垂直方向线
load ([g_src 'LaserHeightH.mat'], 'LaserHeightH');%laser2 是水平方向线

[position_2D1, position_3D1] = getLaserHeight2D(LaserHeightV, rect);
[position_2D2, position_3D2] = getLaserHeight2D(LaserHeightH, rect);
position_2D = [position_2D1; position_2D2]';
position_3D = [position_3D1; position_3D2]';

% 显示激光高度
laserH = displayLaserH(position_2D,position_3D(3,:),rectImg);

%% -- step 2 -- 光度立体
% 计算目标法向（梯度）、初始高度
[p,q,mask_ps] = getGradients(g_src, rect);

Height_original = poisson_solver_function_neumann(p,q);
figure;mesh(p); figure;mesh(q);

Height_original(~mask_ps) = NaN;
figure;mesh(Height_original);
set(gcf,'unit','centimeters','position',[20 5 11 6]);
set(gca,'Position',[.1 .1 .88 .88]); 
xlabel('X(pixel)'),ylabel('Y(pixel)'),zlabel('Z(pixel)'); axis equal;

%% 局部物体的激光三维有效点筛选。如贝壳，则使用此功能。
if control_part == 1
    Height_original = wls2015integration(p,q,mask_ps); Height_original(~mask_ps) = NaN;
    index_2D = position_2D(3,:);
    index_mask = mask_ps(index_2D);
    position_2D = position_2D(:,index_mask);
    position_3D = position_3D(:,index_mask);
end

% realMaskH = imcrop(LaserHeightH.mask, rect).*mask_ps;
% save([g_src 'result\realMaskH.mat'], 'realMaskH');

%%  -- step 3 --若激光有效点数量不足100，可进行梯度偏差自校正
if (length(position_3D(1,:)) < 100) && (control_part == 0)
    disp('step 3 --梯度偏差自校正');
    %%% 方法1：原始梯度自校正
%     [p_c, p_para] = correctGradient(p,0);
%     [q_c, q_para] = correctGradient(q,1);   
%     Height_c = poisson_solver_function_neumann(p_c,q_c);

    %%% 方法2：原始高度自校正
    Height_c = correctHeight_original(Height_original);  % 拟合 H = a x^2 + b y^2 + c x + d y + e;
    [p_c, q_c] = gradient(Height_c);
    
    figure; mesh(Height_c); title('Self-corrected Pixel Height');
    p = p_c; q = q_c; Height_original = Height_c;
end

%% -- step 4 -- 光度立体矫正 -1- 统一尺度 
% 寻找一个转换点,激光高度，将世界坐标系深度转换为像素高度
disp('step 4.1 --统一尺度');
[namda, pos_key, index_key] = findnamda(K, position_2D, position_3D);
% laserH = mean(position_3D(3,:))

% 激光高度从真实值到像素值
laserHeight_pixel = K(1,1) - position_3D(3,:)/namda; %坐标系转换

% 显示激光高度 %% 验证转换是否正确
displayLaserH(position_2D,-laserHeight_pixel,rectImg);

%% 对比方法
%% 原始方法
% Height_c = Height_original;
%% 方法七：梯度自校正优化方法（对比方法）
% Height_c = correctHeight_OGSC(Height_original, p, q, position_2D, laserHeight_pixel, pos_key, mask_ps);
% figure;mesh(Height_c);title('method 7');

%% 方法六：TPS（对比方法）
% Height_c = correctHeight_TPS(position_2D, pos_key, Height_original, laserHeight_pixel,mask_ps);
% figure;mesh(Height_c);title('method 6');

%% 方法四，整体拟合（对比方法） f测 = namda*f真 + a x^2 + b y^2 + c x + d y + e （失败）
% Height_Method1 = correctHeight_Method1(position_2D,pos_key,Height_original,laserHeight_pixel);
% figure;mesh(Height_Method1);title('Height-Method1');

%% -- step 4 -- 光度立体矫正 -2- 光度立体矫正
% 方法三: 该方法必须是十字激光线，有一定凹凸的物体。单激光线无法拟合误差分布，接近平面的物体，尺度估计不准。
disp('step 4.2 --光度立体矫正'); % 浅浮雕的尺度矫正
% diff = 10; scale_new = 2; diff_scale = 10; count = 0; scale = 1;
% disp('count,scale,diff,scale_new,diff_scale');
% disp([count,scale,diff,scale_new,diff_scale]);
% while (diff > 0.9) && (abs(scale_new - 1) > 0.001)
%     count = count + 1;
%     if count == 3
%         break;
%     end
%     scale_old = scale_new;
%     if diff_scale<0.001 || count > 10
%         break;
%     end
%     Height_c = correctHeight_Method2(scale, position_2D, pos_key, Height_original, laserHeight_pixel,mask_ps);
% %     Height_test = Height_c;
% %     Height_test(index) = laserHeight_pixel;
% %     figure;mesh(Height_test);title('Height_test');
%     
% %     scale_new = psHeight_pixel(:)\laserHeight_pixel(:);
%     scale_new = findScale(position_2D,laserHeight_pixel,Height_c);
%     
%     index = position_2D(3,:)';
%     psHeight_pixel = Height_c(index);
%     diff = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)));  
%     diff_scale = abs(scale_new - scale_old);
%     scale = scale * scale_new;
%     %diff 像素的误差
%     display([count,scale,diff,scale_new,diff_scale]);
% end

% 方法三-modify:实验发现，无需迭代矫正，一次矫正效果就较好。
scale = 1; 
Height_c = correctHeight_Method2(scale, position_2D, pos_key, Height_original, laserHeight_pixel,mask_ps);

scale_new = findScale(position_2D,laserHeight_pixel,Height_c);

index = position_2D(3,:)';
psHeight_pixel = Height_c(index);
diff = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)));
diff_scale = scale_new - scale;

disp('count,scale,diff,scale_new,diff_scale');
display([1,scale,diff,scale_new,diff_scale]);

scale = scale_new * scale;
Height_c = correctHeight_Method2(scale, position_2D, pos_key, Height_original, laserHeight_pixel,mask_ps);

figure;mesh(Height_c);title('method 3 -- pixel height');

%%% 验证转换是否正确
% 显示激光高度
index = position_2D(3,:);
displayLaserH(position_2D,-Height_c(index),rectImg);

%% -- step 4 -- 光度立体矫正 -3- 计算真实尺度
disp('step 4.3 --计算真实尺度');
% Height_c = Height_c - Height_c(pos_key); %保持矫正后的高度，不再调整，因此删去
depth_real = -(Height_c - K(1,1))* namda;%坐标系转换
realPosition = getRealAirPosition(K, depth_real);

%% 记录真实高度数据
PS = -(depth_real - depth_real(pos_key));
% PSRL = PS;
% save([g_src 'result\PSRL.mat'], 'PSRL');

%% 记录真实激光高度
index = position_2D(3,:);
mask = zeros(size(rectImg));
mask(index) = 1;
laserH = mask;
laserH(index) = position_3D(3,:);%
temp = position_3D(3,index_key);
laserH = -(laserH - temp);
laserH(mask == 0) = NaN;
figure;mesh(laserH);title('laserH');colorbar; 
% save([g_src 'result\laserH.mat'], 'laserH');

%% 高度平均误差
index = position_2D(3,:)';
psHeight_real = depth_real(index);
laserHeight_real = position_3D(3,:);
diff = mean(abs(psHeight_real(:) - laserHeight_real(:)));
disp(['real height diff = ' num2str(diff)]);

%% 三维平均误差
index = position_2D(3,:)';
ps_realPosition = position_3D;
ps_realPosition(1,:) = realPosition.x(index);
ps_realPosition(2,:) = realPosition.y(index);
ps_realPosition(3,:) = realPosition.z(index);
pos_temp = (ps_realPosition - position_3D).*(ps_realPosition - position_3D);
diff_temp = sqrt(pos_temp(1,:) + pos_temp(2,:) + pos_temp(3,:));
diff = mean(diff_temp);
disp(['real 3D diff = ' num2str(diff)]);

%% mesh显示 % realPosition.y由正变负后，光照位置不对，不知道camlight如何调节！
figure;mesh(realPosition.x, -realPosition.y, -realPosition.z-min(min(-realPosition.z))); %,'EdgeColor', 'None', 'FaceColor', [1 0.5 0.5]
%%设置三维曲面x轴，y轴，z轴，标题对应内容及三个坐标轴的取值范围%%
set(gcf,'unit','centimeters','position',[3 5 12 6]);
set(gca,'Position',[.1 .1 .88 .88]); 
xlabel('X(mm)'),ylabel('Y(mm)'),zlabel('Z(mm)'); axis equal; camlight(45,-75); colorbar;  %camlight(45,-75);title('mesh三维曲面图');caxis([0,10]);

%% 点云显示
Position = [realPosition.x(:) -realPosition.y(:) -realPosition.z(:)];
ptCloud = pointCloud(Position);
%纹理
image = imread([g_src '0.bmp']); texture = imcrop(image,rect);
% ptCloud.Color = [texture(:) texture(:) texture(:)];
ptCloud.Color = reshape(texture(:),size(texture,1)*size(texture,2),3);
figure; pcshow(ptCloud); 
xlabel('X(mm)'),ylabel('Y(mm)'),zlabel('Z(mm)');  %camlight(30,-60); title('三维点云图');
set(gcf,'unit','centimeters','position',[3 5 11 6])
set(gca,'Position',[.08 .08 .9 .9]);

%%——————————————————————————————%%
function [namda1, pos, index] = findnamda(K, position_2D, position_3D)
%%%
% position_2D 包含x = col，y = row，index，三个维度信息
%%%
if size(position_2D,1) ~= 3
    position_2D = position_2D';
end
if size(position_3D,1) ~= 3
    position_3D = position_3D';
end

p_2D(1,:) = position_2D(1,:) - K(1,3);
p_2D(2,:) = position_2D(2,:) - K(2,3);
sum_2D = sum(p_2D.*p_2D, 1); %对列求和
sort_2D = sort(sum_2D);
index_2D = find(sum_2D == sort_2D(floor(length(sum_2D)/8)+1)); % 接近中心一点
index= index_2D(1);

pos = position_2D(3,index);
col = position_2D(1,index); %x = col
row = position_2D(2,index); %y = row
% pos = sub2ind(size(height), row, col);  %关键函数sub2ind
namda1 = abs(position_3D(1,index)./(col - K(1,3))); %缩放倍数,以x轴做缩放
namda2 = abs(position_3D(2,index)./(row - K(2,3))); %缩放倍数,以y轴做缩放
namda3 = abs(position_3D(3,index)./ K(1,1)); %缩放倍数,以z轴做缩放
end
