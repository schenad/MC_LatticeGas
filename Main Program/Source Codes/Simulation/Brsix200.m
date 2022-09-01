%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

MX=[1:800];                                                      %x-axis
MY=[6^(1/2):6^(1/2):100*6^(1/2)];                                %y-axis

interval=1
finaldata=100
aviname='Brsix.avi';

%read the .txt from the folder
load(['Brp.txt']);

for i=0:interval:finaldata-1;
str=['X',num2str(i)];
eval([(str),'=Brp(200*i+1:200*i+200,:)']);
end

clc;

%change .txt into image and save as .jpg
for i=0:interval:finaldata-1;
str=['X',num2str(i)];    
[m,n]=size(eval(str));                                           %calculate the size of matrix
figure
hold on
markersizeA=round(200*sqrt(100/(m*n)));
markersizeB=round(300*sqrt(100/(m*n)));                            %background

%**************************************************************************

[q,p]=find(eval(str)==0);
pp=p*2+q-2;                                                        %function to calculate relationship
plot(pp,q,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size

[q,p]=find(eval(str)==1);
pp=p*2+q-2;
plot(pp,q,'.','Color',[255/255,255/255,0/255],...         %yellow
                  'markersize',markersizeA)
              

%**************************************************************************

axis([0 600 0 200])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[100,100,1000,577]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

set(gca,'LooseInset',get(gca,'TightInset'))                 %Tight in the page

set (gcf,'PaperPositionMode','auto')                       %for save, otherwise save image will change the size
print('-djpeg',num2str(i))

close all
end

%change .jpg into a moive
aviobj = avifile(aviname);
aviobj.Quality = 100;
aviobj.fps=4;
aviobj.compression='None';

for i=0:interval:finaldata-1;
	frame=im2frame(imread(strcat(num2str(i),'.jpg')));
	aviobj = addframe(aviobj,frame);
end

aviobj = close(aviobj);