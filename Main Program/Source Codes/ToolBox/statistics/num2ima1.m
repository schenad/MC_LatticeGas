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
markersizeA=round(160*sqrt(100/(m*n)));
markersizeB=round(150*sqrt(100/(m*n)));                            %background

%**************************************************************************

[q,p]=find(eval(str)==0);
pp=p*2+q-2;                                                        %function to calculate relationship
ppp=pp+200;
ppp=ppp-400*(ppp>400);
plot(pp,q,'.','Color',[254/255,254/255,254/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,q,'.','Color',[254/255,254/255,254/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)
              
[q,p]=find(eval(str)==1);
pp=p*2+q-2;
ppp=pp+200;
ppp=ppp-400*(ppp>400);
plot(pp,q,'.','Color',[0/255,0/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,q,'.','Color',[0/255,0/255,0/255],...         
                  'markersize',markersizeA)              

              
[q,p]=find(eval(str)==0);
qq=q+100;
pp=p*2+q-2+100;                                                        %function to calculate relationship
ppp=pp+200;
ppp=ppp-400*(ppp>400);
plot(pp,qq,'.','Color',[254/255,254/255,254/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,qq,'.','Color',[254/255,254/255,254/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)
              
[q,p]=find(eval(str)==1);
qq=q+100;
pp=p*2+q-2+100;
ppp=pp+200;
ppp=ppp-400*(ppp>400);
plot(pp,qq,'.','Color',[0/255,0/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,qq,'.','Color',[0/255,0/255,0/255],...         
                  'markersize',markersizeA)              
              
              
              
%**************************************************************************

axis([0 346 0 200])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[100,50,700,700]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

set(gca,'LooseInset',get(gca,'TightInset'));                 %Tight in the page

set (gcf,'PaperPositionMode','auto');                       %for save, otherwise save image will change the size

set (gca,'position',[0,0,1,1])

print('-djpeg',num2str(i+1))

close all
end