%%�����߼�������������߶Ƚ���
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   Modify: 2021.4.29
clear; close all; clc;
addpath(genpath(pwd));  % ��ӵ�ǰ·���µ�������Ŀ¼

%% ��ʼ��
src = '..\1�����������궨\calibration20210425-paper\';
g_src = '..\0ʵ������\20210425-paper\4texture_stone2\'; 
%1cushion;2shell;3texture_stone1;4texture_stone2;5throw_pillow1;6throw_pillow2;7word_stone;8sailboat

load([src 'cameraParams']);
K = (cameraParams.IntrinsicMatrix)';
control_part = 1;  %0Ϊȫ����1Ϊ�ֲ�

image = imread([g_src '0.bmp']); img = image(:,:,1);
ROI = 6; sizeImg = floor(size(img)/ROI);
rect = [sizeImg(2)+1, sizeImg(1)+1, sizeImg(2)*4-1, sizeImg(1)*4-1]; %imcrop(img, rect); ����xy��ü�
rectImg = imcrop(img, rect);
clear img image

% ������Ч�����޸��������
K(1,3) = K(1,3) - sizeImg(2)*1; %x -- col
K(2,3) = K(2,3) - sizeImg(1)*1; %y -- row
clear sizeImg

%% -- step 1 -- ��ȡ�������ʵ��ά����
disp('step 1 -- ���漤���߻�ȡ');
load ([g_src 'LaserHeightV.mat'], 'LaserHeightV');%laser1 �Ǵ�ֱ������
load ([g_src 'LaserHeightH.mat'], 'LaserHeightH');%laser2 ��ˮƽ������

[position_2D1, position_3D1] = getLaserHeight2D(LaserHeightV, rect);
[position_2D2, position_3D2] = getLaserHeight2D(LaserHeightH, rect);
position_2D = [position_2D1; position_2D2]';
position_3D = [position_3D1; position_3D2]';

% ��ʾ����߶�
laserH = displayLaserH(position_2D,position_3D(3,:),rectImg);

%% -- step 2 -- �������
% ����Ŀ�귨���ݶȣ�����ʼ�߶�
[p,q,mask_ps] = getGradients(g_src, rect);

Height_original = poisson_solver_function_neumann(p,q);
figure;mesh(p); figure;mesh(q);

Height_original(~mask_ps) = NaN;
figure;mesh(Height_original);
set(gcf,'unit','centimeters','position',[20 5 11 6]);
set(gca,'Position',[.1 .1 .88 .88]); 
xlabel('X(pixel)'),ylabel('Y(pixel)'),zlabel('Z(pixel)'); axis equal;

%% �ֲ�����ļ�����ά��Ч��ɸѡ���籴�ǣ���ʹ�ô˹��ܡ�
if control_part == 1
    Height_original = wls2015integration(p,q,mask_ps); Height_original(~mask_ps) = NaN;
    index_2D = position_2D(3,:);
    index_mask = mask_ps(index_2D);
    position_2D = position_2D(:,index_mask);
    position_3D = position_3D(:,index_mask);
end

% realMaskH = imcrop(LaserHeightH.mask, rect).*mask_ps;
% save([g_src 'result\realMaskH.mat'], 'realMaskH');

%%  -- step 3 --��������Ч����������100���ɽ����ݶ�ƫ����У��
if (length(position_3D(1,:)) < 100) && (control_part == 0)
    disp('step 3 --�ݶ�ƫ����У��');
    %%% ����1��ԭʼ�ݶ���У��
%     [p_c, p_para] = correctGradient(p,0);
%     [q_c, q_para] = correctGradient(q,1);   
%     Height_c = poisson_solver_function_neumann(p_c,q_c);

    %%% ����2��ԭʼ�߶���У��
    Height_c = correctHeight_original(Height_original);  % ��� H = a x^2 + b y^2 + c x + d y + e;
    [p_c, q_c] = gradient(Height_c);
    
    figure; mesh(Height_c); title('Self-corrected Pixel Height');
    p = p_c; q = q_c; Height_original = Height_c;
end

%% -- step 4 -- ���������� -1- ͳһ�߶� 
% Ѱ��һ��ת����,����߶ȣ�����������ϵ���ת��Ϊ���ظ߶�
disp('step 4.1 --ͳһ�߶�');
[namda, pos_key, index_key] = findnamda(K, position_2D, position_3D);
% laserH = mean(position_3D(3,:))

% ����߶ȴ���ʵֵ������ֵ
laserHeight_pixel = K(1,1) - position_3D(3,:)/namda; %����ϵת��

% ��ʾ����߶� %% ��֤ת���Ƿ���ȷ
displayLaserH(position_2D,-laserHeight_pixel,rectImg);

%% �Աȷ���
%% ԭʼ����
% Height_c = Height_original;
%% �����ߣ��ݶ���У���Ż��������Աȷ�����
% Height_c = correctHeight_OGSC(Height_original, p, q, position_2D, laserHeight_pixel, pos_key, mask_ps);
% figure;mesh(Height_c);title('method 7');

%% ��������TPS���Աȷ�����
% Height_c = correctHeight_TPS(position_2D, pos_key, Height_original, laserHeight_pixel,mask_ps);
% figure;mesh(Height_c);title('method 6');

%% �����ģ�������ϣ��Աȷ����� f�� = namda*f�� + a x^2 + b y^2 + c x + d y + e ��ʧ�ܣ�
% Height_Method1 = correctHeight_Method1(position_2D,pos_key,Height_original,laserHeight_pixel);
% figure;mesh(Height_Method1);title('Height-Method1');

%% -- step 4 -- ���������� -2- ����������
% ������: �÷���������ʮ�ּ����ߣ���һ����͹�����塣���������޷�������ֲ����ӽ�ƽ������壬�߶ȹ��Ʋ�׼��
disp('step 4.2 --����������'); % ǳ����ĳ߶Ƚ���
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
%     %diff ���ص����
%     display([count,scale,diff,scale_new,diff_scale]);
% end

% ������-modify:ʵ�鷢�֣��������������һ�ν���Ч���ͽϺá�
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

%%% ��֤ת���Ƿ���ȷ
% ��ʾ����߶�
index = position_2D(3,:);
displayLaserH(position_2D,-Height_c(index),rectImg);

%% -- step 4 -- ���������� -3- ������ʵ�߶�
disp('step 4.3 --������ʵ�߶�');
% Height_c = Height_c - Height_c(pos_key); %���ֽ�����ĸ߶ȣ����ٵ��������ɾȥ
depth_real = -(Height_c - K(1,1))* namda;%����ϵת��
realPosition = getRealAirPosition(K, depth_real);

%% ��¼��ʵ�߶�����
PS = -(depth_real - depth_real(pos_key));
% PSRL = PS;
% save([g_src 'result\PSRL.mat'], 'PSRL');

%% ��¼��ʵ����߶�
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

%% �߶�ƽ�����
index = position_2D(3,:)';
psHeight_real = depth_real(index);
laserHeight_real = position_3D(3,:);
diff = mean(abs(psHeight_real(:) - laserHeight_real(:)));
disp(['real height diff = ' num2str(diff)]);

%% ��άƽ�����
index = position_2D(3,:)';
ps_realPosition = position_3D;
ps_realPosition(1,:) = realPosition.x(index);
ps_realPosition(2,:) = realPosition.y(index);
ps_realPosition(3,:) = realPosition.z(index);
pos_temp = (ps_realPosition - position_3D).*(ps_realPosition - position_3D);
diff_temp = sqrt(pos_temp(1,:) + pos_temp(2,:) + pos_temp(3,:));
diff = mean(diff_temp);
disp(['real 3D diff = ' num2str(diff)]);

%% mesh��ʾ % realPosition.y�����为�󣬹���λ�ò��ԣ���֪��camlight��ε��ڣ�
figure;mesh(realPosition.x, -realPosition.y, -realPosition.z-min(min(-realPosition.z))); %,'EdgeColor', 'None', 'FaceColor', [1 0.5 0.5]
%%������ά����x�ᣬy�ᣬz�ᣬ�����Ӧ���ݼ������������ȡֵ��Χ%%
set(gcf,'unit','centimeters','position',[3 5 12 6]);
set(gca,'Position',[.1 .1 .88 .88]); 
xlabel('X(mm)'),ylabel('Y(mm)'),zlabel('Z(mm)'); axis equal; camlight(45,-75); colorbar;  %camlight(45,-75);title('mesh��ά����ͼ');caxis([0,10]);

%% ������ʾ
Position = [realPosition.x(:) -realPosition.y(:) -realPosition.z(:)];
ptCloud = pointCloud(Position);
%����
image = imread([g_src '0.bmp']); texture = imcrop(image,rect);
% ptCloud.Color = [texture(:) texture(:) texture(:)];
ptCloud.Color = reshape(texture(:),size(texture,1)*size(texture,2),3);
figure; pcshow(ptCloud); 
xlabel('X(mm)'),ylabel('Y(mm)'),zlabel('Z(mm)');  %camlight(30,-60); title('��ά����ͼ');
set(gcf,'unit','centimeters','position',[3 5 11 6])
set(gca,'Position',[.08 .08 .9 .9]);

%%������������������������������������������������������������%%
function [namda1, pos, index] = findnamda(K, position_2D, position_3D)
%%%
% position_2D ����x = col��y = row��index������ά����Ϣ
%%%
if size(position_2D,1) ~= 3
    position_2D = position_2D';
end
if size(position_3D,1) ~= 3
    position_3D = position_3D';
end

p_2D(1,:) = position_2D(1,:) - K(1,3);
p_2D(2,:) = position_2D(2,:) - K(2,3);
sum_2D = sum(p_2D.*p_2D, 1); %�������
sort_2D = sort(sum_2D);
index_2D = find(sum_2D == sort_2D(floor(length(sum_2D)/8)+1)); % �ӽ�����һ��
index= index_2D(1);

pos = position_2D(3,index);
col = position_2D(1,index); %x = col
row = position_2D(2,index); %y = row
% pos = sub2ind(size(height), row, col);  %�ؼ�����sub2ind
namda1 = abs(position_3D(1,index)./(col - K(1,3))); %���ű���,��x��������
namda2 = abs(position_3D(2,index)./(row - K(2,3))); %���ű���,��y��������
namda3 = abs(position_3D(3,index)./ K(1,1)); %���ű���,��z��������
end
