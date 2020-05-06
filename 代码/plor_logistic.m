clear all;clc;   
hold on  
% set(0,'defaultfigurecolor','w');%���ñ�����ɫΪ��ɫ   
axis([0,4,0,1]); 
grid off   
for a=0:0.001:4; 
    x=[0.1234];  
    u=a;  
    for n=2:150      
        x(n)=u* x(n-1)*(1-x(n-1));          
    end
    for n=100:150              
        plot(a,x(n),'k','markersize',3);                
    end
end
title('\fontsize{10}Logisticӳ�����ͼ');          
xlabel('\fontsize{10}��֧����u'),
ylabel('\fontsize{10}������зֲ�x(n)');          
set(gca,'FontSize',10); % �������ִ�С��ͬʱӰ���������ע��ͼ��������ȡ�