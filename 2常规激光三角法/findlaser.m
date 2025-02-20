function mask = findlaser(img, model)
% model:1--rows--(laser_x), 2--cols-(laser_y)
% img: gray img

%% Ѱ��ÿ��ÿ�е����ֵ
[h, w] = size(img);
maxImg = zeros(size(img));
max_rows = zeros(h,1);
max_cols = zeros(w,1);
if model == 1
%     [Cb, bb] = max(img, [], 2); %ÿ�����
    for i = 1:h
        img_row = img(i,:); %ĳһ��
        max_row = max(img_row(:)); %ĳ�����ֵ
        index = find(img_row == max_row);
        index_max = index(round(length(index)/2));
        maxImg(i,index_max) = max_row;
        
        max_rows(i) = max_row;   % �����޳��쳣ֵ
    end
else
%     [Ca, aa] = max(img, [], 1); %ÿ�����
    for i = 1:w
        img_row = img(:,i); %ĳһ��
        max_row = max(img_row(:)); %ĳ�����ֵ
        index = find(img_row == max_row);
        index_max = index(round(length(index)/2));
        maxImg(index_max,i) = max_row;
        
        max_cols(i) = max_row;   % �����޳��쳣ֵ  
    end
end

%% �޳��쳣ֵ  
% mean_rows = mean(max_rows(:));
% c = std(max_rows(:));
% mask = (max_rows > mean_rows + 3*c) | (max_rows < mean_rows - 3*c);
% mask = (max_rows < mean_rows);
% Max_errordata = max(max(max_rows(mask)));

mask = (maxImg > 100);
% figure; imshow(mask);