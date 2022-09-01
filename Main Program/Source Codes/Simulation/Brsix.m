%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

MX=[1:400];                                                      %x-axis
MY=[3^(1/2):3^(1/2):100*3^(1/2)];                                %y-axis

interval=1
finaldata=100
aviname='Brsix.avi';

%read the .txt from the folder
load(['Br.txt']);

for i=0:interval:finaldata-1;
str=['X',num2str(i)];
eval([(str),'=Br(100*i+1:100*i+100,:)']);
end

clc;

%change .txt into image and save as .jpg
for i=0:interval:finaldata-1;

str=['X',num2str(i)];    
    
[m,n]=size(eval(str));                                           %calculate the size of matrix
figure
hold on
markersizeA=round(140*sqrt(100/(m*n)));
markersizeB=round(150*sqrt(100/(m*n)));                            %background

%**************************************************************************

[q,p]=find(eval(str)==0);
pp=p*2+q-2;                                                        %function to calculate relationship
plot(MX(pp),MY(q),'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size

[q,p]=find(eval(str)==1);
pp=p*2+q-2;
plot(MX(pp),MY(q),'.','Color',[255/255,255/255,0/255],...         %yellow
                  'markersize',markersizeA)
              
[q,p]=find(eval(str)==2);
pp=p*2+q-2;
plot(MX(pp),MY(q),'.','Color',[0/255,0/255,205/255],...           %dark blue
                  'markersize',markersizeA)
              
[q,p]=find(eval(str)==3);
pp=p*2+q-2;
plot(MX(pp),MY(q),'.','Color',[255/255,140/255,0/255],...         %orange
                  'markersize',markersizeA)
              
[q,p]=find(eval(str)==4);
pp=p*2+q-2;
plot(MX(pp),MY(q),'.','Color',[32/255,178/255,170/255],...        %light blue
                  'markersize',markersizeA)
              
[q,p]=find(eval(str)==5);
pp=p*2+q-2;
plot(MX(pp),MY(q),'.','Color',[0/255,0/255,205/255],...           %dark blue
                  'markersize',markersizeA)              

%**************************************************************************

axis([0 310 0 180])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[300,200,750,450]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

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












