%���㼤��ͶӰ�����ά����
% --ͼ���ǽ�����Ķ�άͼ��
% --writed by FanHao, 2019-11-16
% --reviseted by Fanhao, 2020-5-28
clear; close all; clc;
src = '..\1�����������궨\calibration20210425-paper\';
src2 = '..\0ʵ������\20210425-paper\1cushion\'; %1cushion;2shell;3texture_stone1;4texture_stone2;5throw_pillow1;6throw_pillow2;7word_stone;8sailboat
load([src 'cameraParams']);
load([src 'planeParams_horizontal']);
load([src 'planeParams_vertical']); %��pa*x + pb*y + pc*z + pd= 0��

model_h = 2; model_v = 1; % model:1-(Ѱ��ÿ�����), 2-(Ѱ��ÿ�����)
originalImage_v = imread([src2 'laserV.bmp']); %��ֱ
originalImage_h = imread([src2 'laserH.bmp']);

% ˮƽ������
%% 1: ͼ��Ԥ����ͼ�����
undistortedImage_h = undistortImage(originalImage_h, cameraParams);
figure;imshow(undistortedImage_h); title('undistortedImage_h');

%% 2����ȡ�����--��Խ������ͼ�� ����Ҫ���ݼ�����ɫ�޸�findlaser������
%���ڶ��ּ�������ȡ����������Ӧ�û���ѡ������ʵķ�����
mask_h = findlaser(undistortedImage_h(:,:,2), model_h); %��ɫ����
figure;imshow(mask_h); title('mask_h');

%% 3����ȡ������Ӧ����ά���� �� ��ʾ
[row_lasers, col_lasers] = find(mask_h>0); %�õ���������
points_lasers = [col_lasers(:),row_lasers(:)]; %���е�����ת��,2��
laserLocation_h = getLaserH_dictionary_Point(points_lasers, cameraParams, planeParams_horizontal);
LaserHeight_h = zeros(size(mask_h));
LaserHeight_h(mask_h>0) = laserLocation_h.z;
figure; mesh(LaserHeight_h);

% ��ֱ������
%% 1: ͼ��Ԥ����ͼ�����
undistortedImage_v = undistortImage(originalImage_v, cameraParams);
figure;imshow(undistortedImage_v); title('undistortedImage_v');

%% 2����ȡ�����--��Խ������ͼ�� ����Ҫ���ݼ�����ɫ�޸�findlaser������
%���ڶ��ּ�������ȡ����������Ӧ�û���ѡ������ʵķ�����
mask_v = findlaser(undistortedImage_v(:,:,2), model_v);
figure;imshow(mask_v); title('mask_v');

%% 3����ȡ������Ӧ����ά���� �� ��ʾ
[row_lasers, col_lasers] = find(mask_v>0); %�õ���������
points_lasers = [col_lasers(:),row_lasers(:)]; %���е�����ת��,2��
laserLocation_v = getLaserH_dictionary_Point(points_lasers, cameraParams, planeParams_vertical);
LaserHeight_v = zeros(size(mask_v));
LaserHeight_v(mask_v>0) = laserLocation_v.z;
figure; mesh(LaserHeight_v);

%% 4����ʾ�뱣��
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
save([src2 'LaserHeightV.mat'],'LaserHeightV'); %laser1 �Ǵ�ֱ������
save([src2 'LaserHeightH.mat'],'LaserHeightH'); %laser2 ��ˮƽ������

% ���������ؽ������ʾ
LaserHeight_h(mask_h == 0) = NaN;
LaserHeight_v(mask_v == 0) = NaN;

LaserHeight_h(mask_v) = LaserHeight_v(mask_v);
figure;  mesh(-LaserHeight_h,'LineWidth',3);
