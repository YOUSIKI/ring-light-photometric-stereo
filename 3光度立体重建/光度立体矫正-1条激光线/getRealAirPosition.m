function realPosition = getRealAirPosition(K, depth_real)
% ͼ���һ������
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   
[h, w] = size(depth_real); %�궨��ͼ��ߴ�
[X,Y]= meshgrid(1:w, 1:h);
x = (X - K(1,3))/K(1,1); %alpha_c = 0
y = (Y - K(2,3))/K(2,2);

realPosition.x = x.*depth_real;
realPosition.y = y.*depth_real;
realPosition.z = depth_real;

