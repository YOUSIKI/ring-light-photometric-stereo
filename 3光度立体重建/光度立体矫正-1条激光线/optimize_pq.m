function [result_p, result_q] = optimize_pq(p, q, para)
%%%
% ÄâºÏ y = a x.^2 + b y.^2 + c x y + d x + e y + f;
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   
%%%
%% p,q ÓÅ»¯
[p_y, p_x] = size(p);
[q_y, q_x] = size(q);

x = 1 : 1 : p_x;
y = 1 : 1 : q_y;

[X,Y]= meshgrid(x,y);

P_1 = 2*para(1)* X + para(3)*Y + para(4);
Q_1 = 2*para(2)* Y + para(3)*X + para(5);

% figure; mesh(p); figure;mesh(P_1);
% figure; mesh(q); figure;mesh(Q_1);

result_p = p - P_1;
result_q = q - Q_1;