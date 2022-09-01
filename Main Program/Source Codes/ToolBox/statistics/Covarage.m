%
clc
clear
close all

%initial: read data and rough process the image
file_name{1}='d83.bmp';
file_name{2}='d88.bmp';
file_name{3}='d93.bmp';
file_name{4}='d94.bmp';
% file_name{5}='d83.bmp';
% file_name{6}='d88.bmp';
% file_name{7}='d93.bmp';
% file_name{8}='d94.bmp';
% file_name{9}='d83.bmp';
% file_name{10}='d88.bmp';

for fi=1:4
    I=imread(file_name{fi});
    thresh=graythresh(I);
    I1=~im2bw(I,thresh);
    p_b(fi)=length(find(I1==1));
    p_t(fi)=size(I1,1)*size(I1,2);
    Coverage(fi)=length(find(I1==1))/size(I1,1)/size(I1,2);
end


T_coverage=mean(Coverage);
Tpoint=size(I1,1)*size(I1,2);
Pt_b=mean(p_b);





