%LR pair correlation function
function [y1,distri,x]=pair_function_dis(data,pix,xy)
%y1: pair 
% distri:potential
% yfit: fitted potential
%xy: size of image


T=273;
Kb=8.6e-5;
% data=xlsread('LR.xls');
ind=1;
for i=1:length(data)
    for j=i:length(data),
        dis(ind)=norm(data(j,:)-data(i,:))*xy/pix;
        ind=ind+1;
    end
end
spa=0.1;
x=0.3:spa:20;
x=x';
y=zeros(length(x),1);
for L=1:length(x)
    a=find(dis>=x(L)&dis<x(L)+spa);
    y(L)=length(a);
end

figure
bar(x,y);hold on
idea1=@(a,x)a*(x.*(1-x.*(4*xy-4*x+pi*x)/(pi*1e4)));
times=20;
times1=nlinfit(x,y,idea1,times);
% times1=5;
idea2=idea1(times1,x);
plot(x,idea2,'r-.');

% 
% figure
y1=y./idea2;
% stem(x,y1,'fill')
% 
% figure
distri=-Kb*T*log(y1);
% stem(x,distri);hold on
% title('pair function')


