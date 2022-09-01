%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Make Image from "*.txt" files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;

siz=100;
        
MX=[1:siz*4];                                                      %x-axis
MY=[3^(1/2):3^(1/2):siz*3^(1/2)];                                %y-axis

interval=1;
finaldata=2;

%read the .txt from the folder
load(['Br.txt']);

for i=0:interval:finaldata-1;
str=['X',num2str(i)];
eval([(str),'=Br(siz*i+1:siz*i+siz,:)']);
end

clc;

%change .txt into image and save as .jpg
for i=0:interval:finaldata-1;
str=['X',num2str(i)];    
[m,n]=size(eval(str));                                           %calculate the size of matrix
figure
hold on
markersizeA=round(siz*sqrt(120/(m*n)));
markersizeB=round(siz*sqrt(100/(m*n)));                            %background

%**************************************************************************

[q,p]=find(eval(str)==0);
pp=p*2+q-2;                                                        %function to calculate relationship
ppp=pp+siz*2;
ppp=ppp-siz*4*(ppp>siz*4);
plot(pp,q,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,q,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)
              
                    
[q,p]=find(eval(str)==0);
qq=q+siz;
pp=p*2+q-2+siz;                                                        %function to calculate relationship
ppp=pp+siz*2;
ppp=ppp-siz*4*(ppp>siz*4);
plot(pp,qq,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)                        %define the dot size
plot(ppp,qq,'.','Color',[105/255,105/255,105/255],...        %'.' means dot, 'color' define gray color
                  'markersize',markersizeB)

              
[q,p]=find(eval(str)==1);
pp=p*2+q-2;
ppp=pp+siz*2;
ppp=ppp-siz*4*(ppp>siz*4);
plot(pp,q,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,q,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)  
              
              
[q,p]=find(eval(str)==1);
qq=q+siz;
pp=p*2+q-2+siz;
ppp=pp+siz*2;
ppp=ppp-siz*4*(ppp>siz*4);
plot(pp,qq,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)
plot(ppp,qq,'.','Color',[255/255,255/255,0/255],...         
                  'markersize',markersizeA)              
              
              
              
%**************************************************************************

axis([0 347*siz/100 0 siz*2])                                       %make sure the white are almost same

axis off;

set (gcf,'Position',[150,50,siz*5.5,siz*5.5]);                    %300,100 is the ordination; 1000,600 is the height and width for window   

set(gca,'LooseInset',get(gca,'TightInset'));                 %Tight in the page

set (gcf,'PaperPositionMode','auto');                       %for save, otherwise save image will change the size

set (gca,'position',[0,0,1,1]);

print('-dpng','-r800',num2str(i+1));

close all;
end