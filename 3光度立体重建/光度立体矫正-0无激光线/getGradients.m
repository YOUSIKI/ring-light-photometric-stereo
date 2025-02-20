function [p,q,mask] = getGradients(g_src, rect)
%PS�ؽ�
%writen by FanHao
%date 2016.8.16

% clear all;close all;clc;
% addpath(genpath(pwd));  % ��ӵ�ǰ·���µ�������Ŀ¼
g_pic_num = 6;
shadowThresh = 0.05;

%% ����ͼƬ���� �� ���շ���
% ͼ������
I = ones(rect(4)+1, rect(3)+1, g_pic_num);
% ���սǶ�
L = zeros(3, g_pic_num);
Slant = 45; %�ӽǷ�����սǶ�
Slant_sin = sin(Slant/180*pi);
Slant_cos = cos(Slant/180*pi);

temp = 0;
for  i = 1:g_pic_num
    %ͼ������
    numberI = (i-1)*360/g_pic_num;
    image = imread([g_src num2str(numberI) '.bmp']);
    image = image(:,:,1);
    image = im2double(image);
    image = imcrop(image,rect); %�ü�Ŀ������
%     imshow(image);
    temp = temp + 1;
    I(:,:, temp) = image /prctile(image(:), 99); %%�о����ǿ���ų��߹�����˼
    %��������--��ʵ�������ϵ
    L(1,temp) = Slant_sin * cos(numberI/180*pi);
    L(2,temp) = Slant_sin * sin(numberI/180*pi);
    L(3,temp) = Slant_cos;
end
clear name Slant Slant_sin Slant_cos i

%% ȥ��Ӱ
% Create a shadow mask.
shadow_mask = (I > shadowThresh);
se = strel('disk', 3);
for i = 1:temp
    % Erode the shadow map to handle de-mosaiking artifact near shadow boundary.
    shadow_mask(:,:,i) = imerode(shadow_mask(:,:,i), se);
end

%% ���㷨����
[rho, n] = PhotometricStereo(I, shadow_mask, L);

%% Visualize the normal map. 
figure; imshow(n);
% figure; imshow(rho); axis xy;
mask = ~isnan(rho); mask = Bigarea(mask);
figure; imshow(mask);

%% Estimate depth map from the normal vectors.
fprintf('Estimating depth map from normal vectors...\n');
p = -n(:,:,1) ./ n(:,:,3);
q = -n(:,:,2) ./ n(:,:,3);
p(isnan(p)) = 0; %�ж������Ԫ���Ƿ���NaN�� NaN �� Not a Number ����д��
q(isnan(q)) = 0;