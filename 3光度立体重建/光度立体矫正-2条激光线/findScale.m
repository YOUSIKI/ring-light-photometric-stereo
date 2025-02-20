function scale = findScale(position_2D,laserHeight_pixel,Height_c)
%
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%  
%% �������������
cols = position_2D(1,:); %x = col
rows = position_2D(2,:); %y = row
index = position_2D(3,:)'; %��Ч����߶ȵ�����
psHeight_pixel = Height_c(index);

x = cols(:);
y = rows(:);
% AA = [psHeight_pixel(:) x, y]; %�г����Ч������
AA = [psHeight_pixel(:) x, y, ones(length(x),1)];
% para_ps = AA \ psHeight_pixel(:);  %�ֿ���ϣ�Ч������
para_laser = AA \ laserHeight_pixel(:);

% ps = psHeight_pixel - AA*para_laser(:);
% laser = laserHeight_pixel - AA*para_laser(:);
% 
% ps = ps - ps(floor(length(laser)/2));
% laser = laser - laser(floor(length(laser)/2));
% scale = ps(:)\laser(:);
scale = para_laser(1);