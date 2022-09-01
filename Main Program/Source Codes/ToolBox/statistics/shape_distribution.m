function shape_island=shape_distribution(I1,Ilabel,size_nm)
%initial: read data and rough process the image

figure
imagesc(I1);title('input figure');

%
% %step 1:shape

    stat = regionprops(Ilabel,'area');
    num_islands=max(max(Ilabel));
    area1=zeros(length(num_islands),1);
    for x = 1: numel(stat)
        area1(x)=[stat(x).Area(1)];
    end

    STATS=regionprops(Ilabel,'Perimeter');
    perimeter1=zeros(length(num_islands),1);
    for x = 1: numel(STATS)
        perimeter1(x)=[STATS(x).Perimeter(1)];
    end
   
    shape_island=(perimeter1.^2)./(area1-perimeter1./2)/4/pi;