function [position_2D, position_3D] = getLaserHeight2D(LaserHeight, rect, mask)
%%%
%position_2D -- x,y,index
%LaserHieht -- x,y,z,mask
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   
%%%
real_x = imcrop(LaserHeight.x, rect);
real_y = imcrop(LaserHeight.y, rect);
real_z = imcrop(LaserHeight.z, rect);

if (~exist('mask','var') || isempty(mask)) 
    mask=ones(size(real_x));
end

realMask = imcrop(LaserHeight.mask, rect).*mask;

[indexY, indexX] = find(realMask == 1);
index = find(realMask == 1);
position_2D = [indexX indexY index];
position_3D = [real_x((realMask == 1)) real_y((realMask == 1)) real_z((realMask == 1))];

% %% 验证转换是否正确
% Height = zeros(size(real_x));
% cols = position_2D(:,1); %x = col
% rows = position_2D(:,2); %y = row
% Height(sub2ind(size(real_x), rows, cols)) = position_3D(:,3);  %关键函数sub2ind
% 
% figure;mesh(real_z);title('real_z');
% figure;mesh(Height);title('Height');
% figure;

