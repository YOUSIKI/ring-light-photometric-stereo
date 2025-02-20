function Height_c = correctHeight_TPS(position_2D, pos_key, Height_original, laserHeight_pixel, mask)

% 1 根据参考点，在像素尺度上做矫正
% 2 采用TPS（薄板样条插值）结合稀疏点将光度立体的结果进行修正
% parameters：
% position_2D：2D position in image of points
% position_3D：3D location related to 2D position of the points
% K: camera intrinsic parameter
% height_ps: height constructed by photometric stereo in pixel size
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   Reference: Depth from gradient fields and control points: bias correction in photometric stereo
%%%%%

%% --1--缩放是近似放大photometric stereo的三维成像，因为是透视投影，缩放存在不可避免的误差。
if size(position_2D,1) ~= 3
    position_2D = position_2D';
end

height_ps = Height_original - Height_original(pos_key);

index = position_2D(3,:)';
key_ps = height_ps(index);  

diff = key_ps(:) - laserHeight_pixel(:);

%% 数据太多，影响TPS计算速度，因此随机选择100个点
len = length(diff);
if len > 100
    rand = randperm(len);
    index_rand = rand(1:100);
end
position_2D = position_2D(:,index_rand);
diff = diff(index_rand);

%% TPS
% st = tpaps(position_2D, diff',1);                %计算形变系数
st = tpaps(position_2D(1:2,:), diff');                %计算形变系数
%  fnplt(st), hold on

[m, n] = size(height_ps);
[x_img, y_img] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系
xy_img = [x_img(:)'; y_img(:)'];
vals = stval(st, xy_img);                         %计算全图的形变对应结果
diff_real = reshape(vals(:), m, n);
Height_c = height_ps - diff_real;
% Height_c = Height_c.*mask;

psHeight_pixel = Height_c(index);% - Height_c(pos_key);
diff_TPS = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)));
display(diff_TPS);