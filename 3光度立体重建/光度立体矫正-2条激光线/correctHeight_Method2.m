function Height = correctHeight_Method2(scale, position_2D, pos_key, Height_original, laserHeight_pixel, mask_ps)
%%%
% 拟合 y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;
%
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%  
%%%
% scale = 1.4;

cols = position_2D(1,:); %x = col
rows = position_2D(2,:); %y = row
index = sub2ind(size(Height_original), rows, cols);
psHeight_pixel = Height_original(index) - Height_original(pos_key);

psHeight_pixel = scale * psHeight_pixel(:) - laserHeight_pixel(:);
x = cols(:);
y = rows(:);
AA = [x.^2, y.^2, x, y, ones(length(x),1)];
% AA = [ x.^2, y.^2, x.*y, x, y, ones(length(x),1)]; %效果不好
%方法1和方法2 结果相同
%% 方法1：矩阵除法
para = AA \ psHeight_pixel;

%% 根据参数计算真实结果
[m, n] = size(Height_original);
[xx, yy] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系
fitted_height = [xx(:).^2, yy(:).^2, xx(:), yy(:), ones(m*n,1)]*para(:);
% ht = reshape(fitted_height, m, n);
% figure; mesh(ht);
% fitted_height = fitted_height - temp_height(pos_key);
Height_temp = Height_original - Height_original(pos_key);

Height = scale * Height_temp(:) - fitted_height(:);
% Height = Height - Height(pos_key); %保持矫正后的高度，不再调整，因此删去
Height = reshape(Height, m, n);

%% 显示拟合激光平面
%显示激光点
if scale == 1
    if length(psHeight_pixel(:)) > 500
        index = randi(length(psHeight_pixel(:)),1,500);
    else
        index = 1:length(psHeight_pixel(:));
    end
    figure, plot3(x(index),y(index),psHeight_pixel(index),'o'); axis equal;
    %显示激光平面
    hold on;
    fitted_height = reshape(fitted_height, m, n);
    fitted_height(~mask_ps) = NaN;
    mesh(xx,yy,fitted_height);
    xlabel('X(pixel)'); ylabel('Y(pixel)'); zlabel('Z(pixel)');
    set(gcf,'unit','centimeters','position',[30 5 11 6])
    grid on;
    hold off;
end