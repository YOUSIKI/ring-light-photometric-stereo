function Height_c = correctHeight_original(Height_original)
%%%
% 拟合 H = a x^2 + b y^2 + c x + d y + e;
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%%%
%% 方法：矩阵除法
[m, n] = size(Height_original);
[x, y] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系

AA = [x(:).*x(:), y(:).*y(:), x(:), y(:), ones(n*m,1)];
para = AA \ Height_original(:);

%% 根据参数计算真实结果
Height_fit = AA*para;
% mesh(reshape(temp, m, n));hold on;
Height_c = Height_original(:) - Height_fit(:);
Height_c = reshape(Height_c, m, n);