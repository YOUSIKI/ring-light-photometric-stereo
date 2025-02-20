function A_BW=Bigarea(A)

%   2017.6.1
%	Author : Hao Fan - Ocean University of China 
%	fanhao@ouc.edu.cn

L = bwlabel(A);
stats = regionprops(L);%?
Ar = cat(1,stats.Area);
ind = find(Ar==max(Ar));
A_BW = zeros(size(A));
A_BW(L==ind) = 1;
A_BW = A_BW > 0;
