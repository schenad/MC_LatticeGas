%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

MX=[1:400];                                                      %x-axis
MY=[3^(1/2):3^(1/2):100*3^(1/2)];                                %y-axis

interval=1;
finaldata=10;

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
plot(pp,q,'.','Color',[254/255,254/255,254/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size

[q,p]=find(eval(str)==1);
pp=p*2+q-2;
plot(pp,q,'.','Color',[0/255,0/255,0/255],...         %yellow
                  'markersize',markersizeA)
              

%**************************************************************************

axis([0 300 0 100])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[300,100,500,289]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

set(gca,'LooseInset',get(gca,'TightInset'))                 %Tight in the page

set (gcf,'PaperPositionMode','auto')                       %for save, otherwise save image will change the size
print('-djpeg',num2str(i+1))

close all
end