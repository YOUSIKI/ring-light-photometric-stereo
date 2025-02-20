function laserH = displayLaserH(position_2D,position_height,rectImg)
%%%% ��ʾ����߶�
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   
index = position_2D(3,:);
mask = zeros(size(rectImg));
mask(index) = 1;
laserH = mask;
% position_height = position_3D(3,:);
position_height = max(position_height) - position_height;
laserH(index) = position_height;%
laserH(mask == 0) = NaN;
figure;mesh(laserH);
title('laserH');
% figure;imshow(mask);

% %% ������ʾ
% realPosition.x = position_3D(1,:);
% realPosition.y = position_3D(2,:);
% realPosition.z = position_height;
% 
% Position = [realPosition.x(:) -realPosition.y(:) realPosition.z(:)];
% ptCloud = pointCloud(Position);
% figure;pcshow(ptCloud);
% xlabel('X(mm)'),ylabel('Y(mm)'),zlabel('Z(mm)'); 
% set(gcf,'unit','centimeters','position',[3 5 10 6])
% set(gca,'Position',[.08 .08 .9 .9]);