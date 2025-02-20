%% ���㼤��ƽ��;vertical;��ֱ
%%%%%%%
clear;close all;clc;
%%���ر궨����
load cameraParams cameraParams;
intrisic = cameraParams.IntrinsicMatrix;
RotationM = cameraParams.RotationMatrices;
TranslationV = cameraParams.TranslationVectors;

%%�涨�궨��������궨ͼ��Χ
point_num = 10;
p = zeros(3,point_num);
Loc = [];

inter = intrisic';
Rc_1 = RotationM(:,:,1)';
Tc_1 = TranslationV(1,:);
Rc_2 = RotationM(:,:,2)';
Tc_2 = TranslationV(2,:);
Rc_3 = RotationM(:,:,3)';
Tc_3 = TranslationV(3,:);
Rc_4 = RotationM(:,:,4)';
Tc_4 = TranslationV(4,:);
Rc_5 = RotationM(:,:,5)';
Tc_5 = TranslationV(5,:);

%%������ȡ��������ά����
p(:,1) = [807;270;1];
temp = laser2Dto3D(p(:,1), inter, Rc_1, Tc_1);
Loc = [Loc temp];

p(:,2) = [789;793;1];
temp = laser2Dto3D(p(:,2), inter, Rc_1, Tc_1);
Loc = [Loc temp];

p(:,3) = [842;347;1];
temp = laser2Dto3D(p(:,3), inter, Rc_2, Tc_2);
Loc = [Loc temp];

p(:,4) = [970;712;1];
temp = laser2Dto3D(p(:,4), inter, Rc_2, Tc_2);
Loc = [Loc temp];

p(:,5) = [987;254;1];
temp = laser2Dto3D(p(:,5), inter, Rc_3, Tc_3);
Loc = [Loc temp];

p(:,6) = [967;641;1];
temp = laser2Dto3D(p(:,6), inter, Rc_3, Tc_3);
Loc = [Loc temp];

p(:,7) = [1090;301;1];
temp = laser2Dto3D(p(:,7), inter, Rc_4, Tc_4);
Loc = [Loc temp];

p(:,8) = [832;800;1];
temp = laser2Dto3D(p(:,8), inter, Rc_4, Tc_4);
Loc = [Loc temp];

p(:,9) = [912;237;1];
temp = laser2Dto3D(p(:,9), inter, Rc_5, Tc_5);
Loc = [Loc temp];

p(:,10) = [887;903;1];
temp = laser2Dto3D(p(:,10), inter, Rc_5, Tc_5);
Loc = [Loc temp];


%%��ϼ�������άƽ�淽�� z = d + a*x + b*y
X = Loc(1,:);
Y = Loc(2,:);
Z = Loc(3,:);

xyz = [ones(point_num,1) X' Y'];
para = regress(Z', xyz);
fprintf('%.6f \n', para);
planeParams_vertical = [para(2), para(3), -1, para(1)];
save planeParams_vertical planeParams_vertical

%% ��ʾ��ϼ���ƽ��
xfit = min(X):1:max(X);  
yfit = min(Y):1:max(Y);
[XFIT,YFIT]= meshgrid(xfit,yfit); 
ZFIT = para(1) + para(2) * XFIT + para(3) * YFIT;

%��ʾ�����
figure, plot3(X,Y,Z,'o'); axis equal;
%��ʾ����ƽ��
hold on; 
mesh(XFIT,YFIT,ZFIT);
xlabel('X'); ylabel('Y'); zlabel('Z')
grid on;
hold off;