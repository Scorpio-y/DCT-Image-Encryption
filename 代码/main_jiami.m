%% ����DCT�任����ͼ��ѹ������ϵͳ
%   �ο���http://blog.csdn.net/ahafg/article/details/48808443
%   @author:����
%   @date:2018.03.22
%-------------------------------------------------------------------------------------------------------%
clear;clc;
Image=imread('../����ͼƬ�ͽ��/lena.png','png');         %��ȡͼ����Ϣ
figure;imshow(Image);title('ԭʼͼƬ');
figure;imhist(Image);title('ԭʼͼ��ֱ��ͼ');
[M,N]=size(Image);                      %��ͼ������и�ֵ��M,N
t=8;    %�ֿ��С

%% ԭʼͼƬR,G,Bͨ����Ϣ��
%Rͨ��
T1=imhist(Image);   %ͳ��ͼ��Rͨ���Ҷ�ֵ��0~255�ķֲ����������T1
S1=sum(T1);     %��������ͼ��Rͨ���ĻҶ�ֵ
xxs1=0;           %ԭʼͼƬRͨ�������
for i=1:256
    pp1=T1(i)/S1;   %ÿ���Ҷ�ֵռ�ȣ���ÿ���Ҷ�ֵ�ĸ���
    if pp1~=0
        xxs1=xxs1-pp1*log2(pp1);
    end
end

%% ����
M1=mod(M,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0
N1=mod(N,t);    %����Ϊ�̶���Կ���Ա����ʱ����ȥ�����ϵ�0
if M1~=0
    Image(M+1:M+t-M1,:)=0;
end
if N1~=0
    Image(:,N+1:N+t-N1)=0;
end
[M,N]=size(Image);  %����������������������N��Ϊ��Կ������ʱ��ͨ��N��֪ԭͼ�ĳ��Ϳ�
SUM=M*N;

%% DCT����
T=dctmtx(t);        %8��DCT�任����
I=im2double(Image);
fun = @(block_struct) T*block_struct.data*T';       %Matlab������blockproc��������blkproc����
J = blockproc(I,[t t],fun);
% J = blkproc(I,[8,8],'P1*x*P2',T,T');  %T'ΪT��ת��

%% ѹ��
% ����ѹ������
rate=5;     %ѹ������1���ڵ�������������Ϊ��Կ
YS=zeros(t,t);  % ѹ������
for i=1:rate
    for j=1:rate+1-i
        YS(i,j)=1;
    end
end
a = sum(YS(:));     %ÿ���ӿ鱣�������ݸ���
fun = @(block_struct) block_struct.data.*YS;        %Matlab������blockproc��������blkproc����
J1 = blockproc(J,[t t],fun);
% J1 = blkproc(J,[8,8],'P1.*x',YS);  %�ӿ��ˣ�����ѹ������ÿ���ӿ鱣��a��Ԫ��
l=(SUM*a)/(t^2);   %���б�����Ԫ�ظ����� 61440
% JJ2  = blkproc(J1,[8,8],'P1*x*P2',T',T);    %�ع�ͼ��
% figure;imshow(JJ2);

%��ȡJ1���Ͻ���ӦԪ�أ���������A
n=1;
A=zeros(1,l);     % Ԥ�����ڴ�
for m=1:M/t
    for mm=1:N/t
        J2=J1(t*(m-1)+1:m*t,t*(mm-1)+1:mm*t);
        for i=1:rate
            for j=1:rate+1-i
                A(n)=J2(i,j);
                n=n+1;
            end
        end
    end
end

%% �õ��û�������Ҫ�ľ���h1
%����Logistic��������p
u0=3.98;     %Logistic�����̣��Զ�Ϊ3.98������Ϊ��Կ
x0=sum(sum(I(1:M/2,:)));     %����ó�Logistic��ֵx0������Ϊ��Կ
x0=(2*x0)/SUM;
x0=floor(x0*10^4)/10^4;     %����4λС��
H1=zeros(1,l+1000);        %Ԥ�����ڴ�
H1(1)=x0;
for i=1:l+999                 %����l+999��ѭ�������õ�l+1000�㣨������ֵ��
    H1(i+1)=u0*H1(i)*(1-H1(i));
end
H1=H1(1001:length(H1));            %ȥ��ǰ1000�㣬��ø��õ������
[~,h1]=sort(H1,'descend');

%% �õ����ż�����Ҫ�ľ���H2
u1=3.892;     %Logistic�����̣��Զ�Ϊ3.892������Ϊ��Կ
x1=sum(sum(I(M/2+1:M,:)));     %����ó�Logistic��ֵx0������Ϊ��Կ
x1=(2*x1)/SUM;
x1=floor(x1*10^4)/10^4;     %����4λС��
H2=zeros(1,l+1000);        %Ԥ�����ڴ�
H2(1)=x1;
for i=1:l+999                 %����L+999��ѭ�������õ�L+1000�㣨������ֵ��
    H2(i+1)=u1*H2(i)*(1-H2(i));
end
H2=H2(1001:length(H2));            %ȥ��ǰ1000�㣬��ø��õ������
%��H2ת����-1��1����
H2(H2>0.5)=1;
H2(H2<=0.5)=-1;

%% �ظ�����count��
count = 100;     %�ظ�����count�Σ��ɸ��ģ�����Ϊ��Կ
B1 = A;
B=zeros(1,l);       %Ԥ�����ڴ�
for i = 1:count
    for j=1:l
        B(j)=B1(h1(j));     %�û�����
    end
    B1=B.*H2;               %���ż���
end
save ../���ܺ�����/jiami.mat B1;
%% ������ ֱ�ӽ�ѹ
% ��һά����ת��Ϊԭ��ÿ���ֿ��Ķ�ά���е����Ͻǲ��֣����ಿ�ֲ���
GG=zeros(M,N);
for i=1:SUM/t^2
    j=1;
    B2=B1(1+a*(i-1):i*a); % ԭ��ÿ���ֿ���ȡ���Ͻ�15��Ԫ��
    G3=zeros(rate);
    for m=1:rate
        for n=1:rate+1-m
            G3(m,n)=B2(j);
            j=j+1;
        end
    end
    %ȷ���˷ֿ�ĺ�������
    xx=floor(i/(N/t))+1;      %�˷ֿ����ڵڼ�����
    yy=mod(i,N/t);            %�˷ֿ�ڼ�����
    if yy==0
        xx=xx-1;
        yy=N/t;
    end
    GG((xx-1)*t+1:(xx-1)*t+rate,(yy-1)*t+1:(yy-1)*t+rate)=G3; % �ϲ����зֿ�
end

fun = @(block_struct) T'*block_struct.data*T;       %Matlab������blockproc��������blkproc����
II = blockproc(GG,[t t],fun);
% II = blkproc(GG,[8,8],'P1*x*P2',T',T);    %�ع�ͼ��
II_1=im2uint8(II);
imwrite(II_1,'../����ͼƬ�ͽ��/������ֱ�ӽ�ѹ���lena.png','png');
figure;imshow(II_1);title('�����ܣ�ֱ�ӽ�ѹ���ͼƬ');
figure;imhist(II_1);title('����ͼ��ֱ��ͼ');

% LL1 = (N*a)/(2*t);
% LL2 = (2*M)/t;
% for i = 1:LL2
%     for j = 1:LL1
%         II(i,j) = B1((i-1)*LL1+j);
%     end
% end
% II_1=II;
% imwrite(II_1,'../����ͼƬ�ͽ��/������ֱ�ӽ�ѹ���lena.png','png');
% figure;imshow(II_1);title('�����ܣ�ֱ�ӽ�ѹ���ͼƬ');
% figure;imhist(II_1);title('����ͼ��ֱ��ͼ');
%% ���ܺ���Ϣ��
%Rͨ��
T2=imhist(II_1);
S2=sum(T2);
xxs2=0;
for i=1:256
    pp2=T2(i)/S2;
    if pp2~=0
        xxs2=xxs2-pp2*log2(pp2);
    end
end

%% ��Ϣ���
disp('���ܳɹ�');
disp('��Կ��');
disp(['��Կ1����0=',num2str(u0),'     ��Կ2��x0=',num2str(x0),'    ��Կ3����1=',num2str(u1),'    ��Կ4��x1=',num2str(x1)]);
disp(['��Կ5��count=',num2str(count),'   ��Կ6��N=',num2str(N),'   ��Կ7��M1=',num2str(M1),'   ��Կ8��N1=',num2str(N1),'   ��Կ9��rate=',num2str(rate)]);
disp(['ԭʼͼƬ��Ϣ��=',num2str(xxs1),'����ͼƬ��Ϣ��=',num2str(xxs2)]);