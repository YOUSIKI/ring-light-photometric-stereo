function laserH_3D = getLaserH_dictionary_Point(points, cameraParams, planeParams)
%根据相机标定信息、激光平面标定信息，得到所有图像激光点对应的三维坐标。
% --输入参数points为已提取的矫正图像的激光点的二维坐标
% --输入参数cameraParams为相机的标定结果
% --输入参数planeParams为标定的激光平面参数 （pa*x + pb*y + pc*z + pd= 0）
% --返回参数laserH_3D为计算得到的相机坐标下的三维坐标

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

%% 计算所有图像激光点对应的三维坐标
disp('--3--compute laser height');
% 图像归一化坐标
X = points(:,1); Y = points(:,2);
x = (X - cc_1)/f_1 - alpha_c*(Y-cc_2)/f_2;
y = (Y - cc_2)/f_2;

% compute laser height
h_z = -pd ./ (pa.*x + pb.*y + pc);
h_x = x .* h_z;
h_y = y .* h_z;
% h_z = pd - h_z; %相机坐标系变换到自定义坐标系

laserH_3D.x = h_x;
laserH_3D.y = h_y;
laserH_3D.z = h_z;
