%计算激光投影点的三维坐标
% --图像都是矫正后的二维图像
% --writed by FanHao, 2019-11-16
% --reviseted by Fanhao, 2020-5-28
clear; close all; clc;
src = '..\1相机及激光面标定\calibration20210425-paper\';
src2 = '..\0实验数据\20210425-paper\1cushion\'; %1cushion;2shell;3texture_stone1;4texture_stone2;5throw_pillow1;6throw_pillow2;7word_stone;8sailboat
load([src 'cameraParams']);
load([src 'planeParams_horizontal']);
load([src 'planeParams_vertical']); %（pa*x + pb*y + pc*z + pd= 0）

model_h = 2; model_v = 1; % model:1-(寻找每行最大), 2-(寻找每列最大)
originalImage_v = imread([src2 'laserV.bmp']); %垂直
originalImage_h = imread([src2 'laserH.bmp']);

% 水平方向处理
%% 1: 图像预处理，图像矫正
undistortedImage_h = undistortImage(originalImage_h, cameraParams);
figure;imshow(undistortedImage_h); title('undistortedImage_h');

%% 2：提取激光点--针对矫正后的图像 （需要根据激光颜色修改findlaser函数）
%存在多种激光线提取方法，根据应用环境选择最合适的方法。
mask_h = findlaser(undistortedImage_h(:,:,2), model_h); %绿色激光
figure;imshow(mask_h); title('mask_h');

%% 3：提取激光点对应的三维坐标 并 显示
[row_lasers, col_lasers] = find(mask_h>0); %得到的是行列
points_lasers = [col_lasers(:),row_lasers(:)]; %行列到坐标转换,2列
laserLocation_h = getLaserH_dictionary_Point(points_lasers, cameraParams, planeParams_horizontal);
LaserHeight_h = zeros(size(mask_h));
LaserHeight_h(mask_h>0) = laserLocation_h.z;
figure; mesh(LaserHeight_h);

% 垂直方向处理
%% 1: 图像预处理，图像矫正
undistortedImage_v = undistortImage(originalImage_v, cameraParams);
figure;imshow(undistortedImage_v); title('undistortedImage_v');

%% 2：提取激光点--针对矫正后的图像 （需要根据激光颜色修改findlaser函数）
%存在多种激光线提取方法，根据应用环境选择最合适的方法。
mask_v = findlaser(undistortedImage_v(:,:,2), model_v);
figure;imshow(mask_v); title('mask_v');

%% 3：提取激光点对应的三维坐标 并 显示
[row_lasers, col_lasers] = find(mask_v>0); %得到的是行列
points_lasers = [col_lasers(:),row_lasers(:)]; %行列到坐标转换,2列
laserLocation_v = getLaserH_dictionary_Point(points_lasers, cameraParams, planeParams_vertical);
LaserHeight_v = zeros(size(mask_v));
LaserHeight_v(mask_v>0) = laserLocation_v.z;
figure; mesh(LaserHeight_v);

%% 4：显示与保存
LaserHeightV.x = zeros(size(mask_v)); LaserHeightV.y = LaserHeightV.x; LaserHeightV.z = LaserHeightV.x;
LaserHeightV.x(mask_v>0) = laserLocation_v.x; 
LaserHeightV.y(mask_v>0) = laserLocation_v.y; 
LaserHeightV.z(mask_v>0) = laserLocation_v.z;
LaserHeightV.mask = mask_v;
LaserHeightH.x = zeros(size(mask_v)); LaserHeightH.y = LaserHeightV.x; LaserHeightH.z = LaserHeightV.x;
LaserHeightH.x(mask_h>0) = laserLocation_h.x; 
LaserHeightH.y(mask_h>0) = laserLocation_h.y; 
LaserHeightH.z(mask_h>0) = laserLocation_h.z;
LaserHeightH.mask = mask_h;
save([src2 'LaserHeightV.mat'],'LaserHeightV'); %laser1 是垂直方向线
save([src2 'LaserHeightH.mat'],'LaserHeightH'); %laser2 是水平方向线

% 激光三角重建结果显示
LaserHeight_h(mask_h == 0) = NaN;
LaserHeight_v(mask_v == 0) = NaN;

LaserHeight_h(mask_v) = LaserHeight_v(mask_v);
figure;  mesh(-LaserHeight_h,'LineWidth',3);
