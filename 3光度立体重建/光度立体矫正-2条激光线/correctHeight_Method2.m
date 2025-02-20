function Height = correctHeight_Method2(scale, position_2D, pos_key, Height_original, laserHeight_pixel, mask_ps)
%%%
% ��� y - namda * H  = a x.^2 + b y.^2 + c x + d y + f;
%
%   2018.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn
%  
%%%
% scale = 1.4;

cols = position_2D(1,:); %x = col
rows = position_2D(2,:); %y = row
index = sub2ind(size(Height_original), rows, cols);
psHeight_pixel = Height_original(index) - Height_original(pos_key);

psHeight_pixel = scale * psHeight_pixel(:) - laserHeight_pixel(:);
x = cols(:);
y = rows(:);
AA = [x.^2, y.^2, x, y, ones(length(x),1)];
% AA = [ x.^2, y.^2, x.*y, x, y, ones(length(x),1)]; %Ч������
%����1�ͷ���2 �����ͬ
%% ����1���������
para = AA \ psHeight_pixel;

%% ���ݲ���������ʵ���
[m, n] = size(Height_original);
[xx, yy] = meshgrid(1:n, 1:m);          %ע���С��У���x��y�Ķ�Ӧ��ϵ
fitted_height = [xx(:).^2, yy(:).^2, xx(:), yy(:), ones(m*n,1)]*para(:);
% ht = reshape(fitted_height, m, n);
% figure; mesh(ht);
% fitted_height = fitted_height - temp_height(pos_key);
Height_temp = Height_original - Height_original(pos_key);

Height = scale * Height_temp(:) - fitted_height(:);
% Height = Height - Height(pos_key); %���ֽ�����ĸ߶ȣ����ٵ��������ɾȥ
Height = reshape(Height, m, n);

%% ��ʾ��ϼ���ƽ��
%��ʾ�����
if scale == 1
    if length(psHeight_pixel(:)) > 500
        index = randi(length(psHeight_pixel(:)),1,500);
    else
        index = 1:length(psHeight_pixel(:));
    end
    figure, plot3(x(index),y(index),psHeight_pixel(index),'o'); axis equal;
    %��ʾ����ƽ��
    hold on;
    fitted_height = reshape(fitted_height, m, n);
    fitted_height(~mask_ps) = NaN;
    mesh(xx,yy,fitted_height);
    xlabel('X(pixel)'); ylabel('Y(pixel)'); zlabel('Z(pixel)');
    set(gcf,'unit','centimeters','position',[30 5 11 6])
    grid on;
    hold off;
end