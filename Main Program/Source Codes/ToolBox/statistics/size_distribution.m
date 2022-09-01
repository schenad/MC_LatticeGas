function size_island=size_distribution(I1,Ilabel,size_nm)
%initial: read data and rough process the image

figure
imagesc(I1);title('input figure');

%
% %step 1:size

    stat = regionprops(Ilabel,'area');
    num_islands=max(max(Ilabel));
    area1=zeros(length(num_islands),1);
    for x = 1: numel(stat)
        area1(x)=[stat(x).Area(1)];
    end
    area1=area1*size_nm^2/size(I1,1)^2;
    size_island=round(area1/(0.44*0.44));