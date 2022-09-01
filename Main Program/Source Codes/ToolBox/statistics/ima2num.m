%
clc
clear
close all

xmax=272;
trials=10;

Br=zeros(trials*xmax,xmax);

for fi=1:trials
    I=imread([num2str(fi),'.bmp']);
    thresh=graythresh(I);
    I1=~im2bw(I,thresh);
    I1=imresize(I1,[4*xmax/sqrt(3) 2*xmax],'bicubic');
    [row,col]=find(I1);

for p=1:size(row);
if mod(row(p),4)==1
   if mod(col(p),2)==1
       x=(col(p)+1)/2+((row(p)+1)/2-1)/2;
       y=fi*xmax-(row(p)+1)/2;
       if x<1
           x=x+xmax;
       elseif x>xmax
           x=x-xmax;
       end
       if y<=fi*xmax && y>(fi-1)*xmax
       Br(y,x)=1;
       end
   end
elseif mod(row(p),4)==3;
   if mod(col(p),2)==0
       x=col(p)/2+((row(p)+1)/2-2)/2;
       y=fi*xmax-(row(p)+1)/2;
       if x<1
           x=x+xmax;
       elseif x>xmax
           x=x-xmax;
       end
       if y<=fi*xmax && y>(fi-1)*xmax
       Br(y,x)=1;
       end
   end
end
end
end

dlmwrite('Br.txt',Br,' ');