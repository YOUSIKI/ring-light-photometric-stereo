function B1 = calibration_parameter(height)
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%
[rows, cols] = size(height);
data = zeros(rows*cols, 3);
count = 0;  % 有效数据的个数
for i = 1:5:cols
    for j = 1:5:rows
        count = count+1;
        data(count, 1) = i;
        data(count, 2) = j;
        data(count, 3) = height(j, i);
    end
end

x = data(1:count,1);
y = data(1:count,2);
Z = data(1:count,3);

mask = (~isnan(Z));
x = x(mask);
y = y(mask);
Z = Z(mask);

X = [ x.^2, y.^2, x.*y, x, y,ones(length(Z),1)];
%方法1和方法2 结果相同
%% 方法1：矩阵除法
B1 = X \ Z;
%% 方法2：回归
[B2,BINT,R,RINT,STATS] = regress(Z,X);       %利用回归函数用6点拟合二元二次函数


B = B1;
% backscatter = (B(1) + B(2)*data(:,1).*data(:,1) + B(3)*data(:, 1).*data(:, 2) +...
%     B(4)*data(:, 2).*data(:, 2) + B(5)*data(:, 1) + B(6)*data(:,2));

% 
% diff = data(:,3) - backscatter;

bs = height;
diff = height;
for i = 1:cols
    for j = 1:rows
        bs_x = i;
        bs_y = j;
        bs(j, i) = B(1) + B(2)*bs_x*bs_x + B(3)*bs_y*bs_y + B(4)*bs_x*bs_y + B(5)*bs_x + B(6)*bs_y;

%         diff(j,i) = height(j,i) - bs(j, i);
    end
end
% figure;mesh(bs);
% figure;mesh(diff);