%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

MX=[1:800];                                                      %x-axis
MY=[6^(1/2):6^(1/2):100*6^(1/2)];                                %y-axis

interval=1;
finaldata=10;

%read the .txt from the folder
load(['Br.txt']);

for i=0:interval:finaldata-1;
str=['X',num2str(i)];
eval([(str),'=Br(200*i+1:200*i+200,:)']);
end

clc;

%change .txt into image and save as .jpg
for i=0:interval:finaldata-1;
str=['X',num2str(i)];    
[m,n]=size(eval(str));                                           %calculate the size of matrix
figure
hold on
markersizeA=round(120*sqrt(100/(m*n)));
markersizeB=round(120*sqrt(100/(m*n)));                            %background

%**************************************************************************

[q,p]=find(eval(str)==1);
pp=p*2+q-2;                                                        %function to calculate relationship
ppp=pp+400;
ppp=ppp-800*(ppp>800);
plot(pp,q,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,q,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)
              
[q,p]=find(eval(str)==0);
pp=p*2+q-2;
ppp=pp+400;
ppp=ppp-800*(ppp>800);
plot(pp,q,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,q,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)              

              
[q,p]=find(eval(str)==1);
qq=q+200;
pp=p*2+q-2+200;                                                        %function to calculate relationship
ppp=pp+400;
ppp=ppp-800*(ppp>800);
plot(pp,qq,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,qq,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)
              
[q,p]=find(eval(str)==0);
qq=q+200;
pp=p*2+q-2+200;
ppp=pp+400;
ppp=ppp-800*(ppp>800);
plot(pp,qq,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,qq,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)              
              
              
              
%**************************************************************************

axis([0 693 0 400])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[100,50,700,700]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

set(gca,'LooseInset',get(gca,'TightInset'));                 %Tight in the page

set (gcf,'PaperPositionMode','auto');                       %for save, otherwise save image will change the size

set (gca,'position',[0,0,1,1]);

print('-djpeg','-r800',num2str(i+1));

close all;
end