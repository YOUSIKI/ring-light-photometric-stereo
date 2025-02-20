function Height_c = correctHeight_TPS(position_2D, pos_key, Height_original, laserHeight_pixel, mask)

% 1 ���ݲο��㣬�����س߶���������
% 2 ����TPS������������ֵ�����ϡ��㽫�������Ľ����������
% parameters��
% position_2D��2D position in image of points
% position_3D��3D location related to 2D position of the points
% K: camera intrinsic parameter
% height_ps: height constructed by photometric stereo in pixel size
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   Reference: Depth from gradient fields and control points: bias correction in photometric stereo
%%%%%

%% --1--�����ǽ��ƷŴ�photometric stereo����ά������Ϊ��͸��ͶӰ�����Ŵ��ڲ��ɱ������
if size(position_2D,1) ~= 3
    position_2D = position_2D';
end

height_ps = Height_original - Height_original(pos_key);

index = position_2D(3,:)';
key_ps = height_ps(index);  

diff = key_ps(:) - laserHeight_pixel(:);

%% ����̫�࣬Ӱ��TPS�����ٶȣ�������ѡ��100����
len = length(diff);
if len > 100
    rand = randperm(len);
    index_rand = rand(1:100);
end
position_2D = position_2D(:,index_rand);
diff = diff(index_rand);

%% TPS
% st = tpaps(position_2D, diff',1);                %�����α�ϵ��
st = tpaps(position_2D(1:2,:), diff');                %�����α�ϵ��
%  fnplt(st), hold on

[m, n] = size(height_ps);
[x_img, y_img] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ
xy_img = [x_img(:)'; y_img(:)'];
vals = stval(st, xy_img);                         %����ȫͼ���α��Ӧ���
diff_real = reshape(vals(:), m, n);
Height_c = height_ps - diff_real;
% Height_c = Height_c.*mask;

psHeight_pixel = Height_c(index);% - Height_c(pos_key);
diff_TPS = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)));
display(diff_TPS);