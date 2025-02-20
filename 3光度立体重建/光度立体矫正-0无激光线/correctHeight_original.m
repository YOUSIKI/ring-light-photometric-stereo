function Height_c = correctHeight_original(Height_original)
%%%
% ��� H = a x^2 + b y^2 + c x + d y + e;
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%%%
%% �������������
[m, n] = size(Height_original);
[x, y] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ

AA = [x(:).*x(:), y(:).*y(:), x(:), y(:), ones(n*m,1)];
para = AA \ Height_original(:);

%% ���ݲ���������ʵ���
Height_fit = AA*para;
% mesh(reshape(temp, m, n));hold on;
Height_c = Height_original(:) - Height_fit(:);
Height_c = reshape(Height_c, m, n);