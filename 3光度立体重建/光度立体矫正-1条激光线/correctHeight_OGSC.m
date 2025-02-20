function Height = correctHeight_OGSC(Height_original, p, q, position_2D, laserHeight_pixel, pos_key, mask)
%%%
% 拟合 y = a x.^2 + b y.^2 + c x y + d x + e y + f;
%
%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%   Reference: Deviation correction method for close-range photometric stereo with nonuniform illumination
%%%
%% 激光高度的有效区域
index = position_2D(3,:)'; %有效激光高度点索引

%% 自校正
scale = 1; p = p*scale; q = q*scale; Height_original = Height_original*scale;
para = calibration_parameter(Height_original);
[p_c, q_c] = optimize_pq(p, q, para);
% Height_c = poisson_solver_function_neumann(p_c,q_c);
Height_c = wls2015integration(p_c, q_c, mask); %带mask的最小二乘权重积分方法
Height_c(~mask) = NaN;
figure; mesh(Height_c);
psHeight_pixel = Height_c(index) - Height_c(pos_key);
diff_A = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)));

%% 自校正优化方法-OGSC
ParaB = zeros(2,2); result = zeros(2,1);
ParaB(1,1) = 2*para(1); ParaB(1,2) = para(3); ParaB(2,2) = 2*para(2); ParaB(2,1) = para(3);
result(1,1) = -para(4); result(2,1) = -para(5);
centerB = ParaB \ result;

% 优化迭代
min_diff = diff_A; min_p = p_c; min_q = q_c; min_interation = 0;
Convergence = 0.5; sscale = 0.1; iteration = 0; % 迭代的次数
while(diff_A > Convergence)
    kkk = 1 + ((-1)^(iteration))* iteration * sscale;
    para(1) = para(1) * kkk;
    para(2) = para(2) * kkk;
    
    para(4) = -(2*para(1) *centerB(1)+ para(3)*centerB(2));
    para(5) = -(2*para(2) *centerB(2)+ para(3)*centerB(1));
    iteration = iteration + 1;
    
    [p_c, q_c] = optimize_pq(p, q, para);
    Height_c = wls2015integration(p_c, q_c, mask); %带mask的最小二乘权重积分方法
%     Height_c = poisson_solver_function_neumann(p_c,q_c);
    psHeight_pixel = Height_c(index) - Height_c(pos_key);
    diff_A = mean(abs(psHeight_pixel(:) - laserHeight_pixel(:)))
    
    if diff_A < min_diff
        min_interation = 1;
        min_diff = diff_A; min_p = p_c; min_q = q_c;
    else
        min_interation = min_interation + 1;
        if min_interation > 5
            break;
        end
    end
end

p_OGSC = min_p; q_OGSC = min_q;
% Height = poisson_solver_function_neumann(p_OGSC, q_OGSC);
Height = wls2015integration(p_OGSC, q_OGSC, mask); %带mask的最小二乘权重积分方法
Height(~mask) = NaN;
diff_OGSC = min_diff;
display(diff_OGSC);