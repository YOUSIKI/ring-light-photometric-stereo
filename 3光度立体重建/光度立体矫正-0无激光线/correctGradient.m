function [gradient_correct, para] = correctGradient(gradient,model)
%%%
% 拟合 y = c x + d y + e;
% model = 0, x轴方向梯度矫正；model = 1，y轴方向梯度矫正。
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%%%
%% 方法：矩阵除法
[m, n] = size(gradient);
[x, y] = meshgrid(1:n, 1:m);          %注意行、列，和x，y的对应关系
if model == 0 %x方向梯度
    AA = [x(:), ones(n*m,1)];
else
    AA = [y(:), ones(n*m,1)];
end
para = AA \ gradient(:);

%% 根据参数计算真实结果
temp = AA*para;
% mesh(reshape(temp, m, n));hold on;
gradient_correct = gradient(:) - temp(:);
gradient_correct = reshape(gradient_correct, m, n);