function [gradient_correct, para] = correctGradient(gradient,model)
%%%
% ��� y = c x + d y + e;
% model = 0, x�᷽���ݶȽ�����model = 1��y�᷽���ݶȽ�����
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%%%
%% �������������
[m, n] = size(gradient);
[x, y] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ
if model == 0 %x�����ݶ�
    AA = [x(:), ones(n*m,1)];
else
    AA = [y(:), ones(n*m,1)];
end
para = AA \ gradient(:);

%% ���ݲ���������ʵ���
temp = AA*para;
% mesh(reshape(temp, m, n));hold on;
gradient_correct = gradient(:) - temp(:);
gradient_correct = reshape(gradient_correct, m, n);