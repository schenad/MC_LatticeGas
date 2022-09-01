%
clc
clear
close all

Num_image=17;
for fi=1:Num_image
    I=imread([num2str(fi),'.png']);
    thresh=graythresh(I);
    I1=~im2bw(I,thresh);
    I1=bwareaopen(I1,3);
    % Ibw = imfill(I1,'holes');
    Ilabel = bwlabel(I1);
    figure
    imagesc(I1);title('input figure');
    % %step 1:Do distance statistics
    stat = regionprops(Ilabel,'centroid');
    num_islands=max(max(Ilabel));
    dis=zeros(length(num_islands),2);
    for x = 1: numel(stat)
        dis(x,:)=[stat(x).Centroid(1),stat(x).Centroid(2)];
    end
    [g_r{fi},distri{fi},x]=pair_function_dis(dis,size(Ilabel,1),76);
    %
   
end

Gt_r=0;
for i=1:Num_image
    Gt_r=Gt_r+g_r{i};
end
Gt_r=Gt_r/Num_image;
figure
plot(x,Gt_r,'rd');
print('-djpeg','G(r)')
close all;