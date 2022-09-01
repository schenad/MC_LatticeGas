trials=10;

for fi=1:trials
    I=imread([num2str(fi),'.bmp']);
    thresh=graythresh(I);
    I1=~im2bw(I,thresh);
    I1=bwareaopen(I1,10);
    % Ibw = imfill(I1,'holes');
    Ilabel = bwlabel(I1);
    size_island=size_distribution(I1,Ilabel,120);
    shape_island=shape_distribution(I1,Ilabel,120);
    
    fid = fopen('size.dat','a+');
    fprintf(fid,'%g\n',size_island);
    fclose(fid);
    
    fid = fopen('shape.dat','a+');
    fprintf(fid,'%g\n',shape_island);
    fclose(fid);
end
close all;