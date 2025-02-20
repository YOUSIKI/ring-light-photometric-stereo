function laserH_3D = getLaserH_dictionary_Point(points, cameraParams, planeParams)
%��������궨��Ϣ������ƽ��궨��Ϣ���õ�����ͼ�񼤹���Ӧ����ά���ꡣ
% --�������pointsΪ����ȡ�Ľ���ͼ��ļ����Ķ�ά����
% --�������cameraParamsΪ����ı궨���
% --�������planeParamsΪ�궨�ļ���ƽ����� ��pa*x + pb*y + pc*z + pd= 0��
% --���ز���laserH_3DΪ����õ�����������µ���ά����

disp('--1--camera calibration info');
% load cameraParams cameraParams;
intrisic = cameraParams.IntrinsicMatrix';
f_1 = intrisic(1,1); f_2 = intrisic(2,2);
cc_1 = intrisic(1,3); cc_2 = intrisic(2,3);
alpha_c = 0; %intrisic(1,2)

disp('--2--laser plane calibration info'); %pa*x + pb*y + pc*z + pd= 0
% load planeParams planeParams;
pa = planeParams(1);
pb = planeParams(2);
pc = planeParams(3);
pd = planeParams(4);

%% ��������ͼ�񼤹���Ӧ����ά����
disp('--3--compute laser height');
% ͼ���һ������
X = points(:,1); Y = points(:,2);
x = (X - cc_1)/f_1 - alpha_c*(Y-cc_2)/f_2;
y = (Y - cc_2)/f_2;

% compute laser height
h_z = -pd ./ (pa.*x + pb.*y + pc);
h_x = x .* h_z;
h_y = y .* h_z;
% h_z = pd - h_z; %�������ϵ�任���Զ�������ϵ

laserH_3D.x = h_x;
laserH_3D.y = h_y;
laserH_3D.z = h_z;
