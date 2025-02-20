function Height = correctHeight_Method1(position_2D, pos_key, Height_original, laserHeight_pixel)
%%%
% ��� y = namda * H + a x.^2 + b y.^2 + c x + d y + 1;
%
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   
%%%

cols = position_2D(1,:); %x = col
rows = position_2D(2,:); %y = row
index = sub2ind(size(Height_original), rows, cols);
psHeight_pixel = Height_original(index) - Height_original(pos_key);

x = cols(:);
y = rows(:);
count = length(x);
AA = [laserHeight_pixel(:) x.^2, y.^2, x, y, ones(count,1)];
% AA = [ones(count,1), x.^2, y.^2, x.*y, x, y];
%����1�ͷ���2 �����ͬ
%% ����1���������
para = AA \ psHeight_pixel(:);

%% ���ݲ���������ʵ���
[m, n] = size(Height_original);
[xx, yy] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ
Height_temp = Height_original - Height_original(pos_key);
figure;mesh(Height_temp);title('Height-temp');
temp_height = [Height_temp(:) xx(:).^2 yy(:).^2 xx(:) yy(:) ones(m*n,1)]*para(:);

Height = temp_height(:) + Height_original(pos_key);
Height = reshape(Height, m, n);