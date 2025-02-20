%% 估算激光平面;horizontal;水平
%%%%%%%
clear;close all;clc;
%%加载标定参数
load cameraParams cameraParams;
intrisic = cameraParams.IntrinsicMatrix;
RotationM = cameraParams.RotationMatrices;
TranslationV = cameraParams.TranslationVectors;

%%规定标定点个数级标定图像范围
point_num = 10;
p = zeros(3,point_num);
Loc = [];

inter = intrisic';
Rc_1 = RotationM(:,:,6)';
Tc_1 = TranslationV(6,:);
Rc_2 = RotationM(:,:,7)';
Tc_2 = TranslationV(7,:);
Rc_3 = RotationM(:,:,8)';
Tc_3 = TranslationV(8,:);
Rc_4 = RotationM(:,:,9)';
Tc_4 = TranslationV(9,:);
Rc_5 = RotationM(:,:,10)';
Tc_5 = TranslationV(10,:);

%%计算提取激光点的三维坐标
p(:,1) = [498;545;1];
temp = laser2Dto3D(p(:,1), inter, Rc_1, Tc_1);
Loc = [Loc temp];

p(:,2) = [1020;552;1];
temp = laser2Dto3D(p(:,2), inter, Rc_1, Tc_1);
Loc = [Loc temp];

p(:,3) = [462;468;1];
temp = laser2Dto3D(p(:,3), inter, Rc_2, Tc_2);
Loc = [Loc temp];

p(:,4) = [1167;391;1];
temp = laser2Dto3D(p(:,4), inter, Rc_2, Tc_2);
Loc = [Loc temp];

p(:,5) = [576;379;1];
temp = laser2Dto3D(p(:,5), inter, Rc_3, Tc_3);
Loc = [Loc temp];

p(:,6) = [1206;479;1];
temp = laser2Dto3D(p(:,6), inter, Rc_3, Tc_3);
Loc = [Loc temp];

p(:,7) = [425;284;1];
temp = laser2Dto3D(p(:,7), inter, Rc_4, Tc_4);
Loc = [Loc temp];

p(:,8) = [1210;479;1];
temp = laser2Dto3D(p(:,8), inter, Rc_4, Tc_4);
Loc = [Loc temp];

p(:,9) = [394;543;1];
temp = laser2Dto3D(p(:,9), inter, Rc_5, Tc_5);
Loc = [Loc temp];

p(:,10) = [1203;555;1];
temp = laser2Dto3D(p(:,10), inter, Rc_5, Tc_5);
Loc = [Loc temp];


%%拟合激光点的三维平面方程 z = d + a*x + b*y
X = Loc(1,:);
Y = Loc(2,:);
Z = Loc(3,:);

xyz = [ones(point_num,1) X' Y'];
para = regress(Z', xyz);
fprintf('%.6f \n', para);
planeParams_horizontal = [para(2), para(3), -1, para(1)];
save planeParams_horizontal planeParams_horizontal

%% 显示拟合激光平面
xfit = min(X):1:max(X);  
yfit = min(Y):1:max(Y);
[XFIT,YFIT]= meshgrid(xfit,yfit); 
ZFIT = para(1) + para(2) * XFIT + para(3) * YFIT;

%显示激光点
figure, plot3(X,Y,Z,'o'); axis equal;
%显示激光平面
hold on; 
mesh(XFIT,YFIT,ZFIT);
xlabel('X'); ylabel('Y'); zlabel('Z')
grid on;
hold off;