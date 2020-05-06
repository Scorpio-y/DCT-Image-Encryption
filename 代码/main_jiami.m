%% 基于DCT变换数字图像压缩加密系统
%   参考：http://blog.csdn.net/ahafg/article/details/48808443
%   @author:沈洋
%   @date:2018.03.22
%-------------------------------------------------------------------------------------------------------%
clear;clc;
Image=imread('../测试图片和结果/lena.png','png');         %读取图像信息
figure;imshow(Image);title('原始图片');
figure;imhist(Image);title('原始图像直方图');
[M,N]=size(Image);                      %将图像的行列赋值给M,N
t=8;    %分块大小

%% 原始图片R,G,B通道信息熵
%R通道
T1=imhist(Image);   %统计图像R通道灰度值从0~255的分布情况，存至T1
S1=sum(T1);     %计算整幅图像R通道的灰度值
xxs1=0;           %原始图片R通道相关性
for i=1:256
    pp1=T1(i)/S1;   %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        xxs1=xxs1-pp1*log2(pp1);
    end
end

%% 补零
M1=mod(M,t);    %可作为固定密钥，以便解码时可以去除补上的0
N1=mod(N,t);    %可作为固定密钥，以便解码时可以去除补上的0
if M1~=0
    Image(M+1:M+t-M1,:)=0;
end
if N1~=0
    Image(:,N+1:N+t-N1)=0;
end
[M,N]=size(Image);  %补零后的行数和列数，其中N作为密钥，解密时可通过N获知原图的长和宽
SUM=M*N;

%% DCT处理
T=dctmtx(t);        %8阶DCT变换矩阵
I=im2double(Image);
fun = @(block_struct) T*block_struct.data*T';       %Matlab建议用blockproc函数代替blkproc函数
J = blockproc(I,[t t],fun);
% J = blkproc(I,[8,8],'P1*x*P2',T,T');  %T'为T的转置

%% 压缩
% 构造压缩矩阵
rate=5;     %压缩矩阵1存在的行列数，可作为密钥
YS=zeros(t,t);  % 压缩矩阵
for i=1:rate
    for j=1:rate+1-i
        YS(i,j)=1;
    end
end
a = sum(YS(:));     %每个子块保留的数据个数
fun = @(block_struct) block_struct.data.*YS;        %Matlab建议用blockproc函数代替blkproc函数
J1 = blockproc(J,[t t],fun);
% J1 = blkproc(J,[8,8],'P1.*x',YS);  %子块点乘，进行压缩，即每个子块保留a个元素
l=(SUM*a)/(t^2);   %所有保留的元素个数和 61440
% JJ2  = blkproc(J1,[8,8],'P1*x*P2',T',T);    %重构图像
% figure;imshow(JJ2);

%提取J1左上角相应元素，生成序列A
n=1;
A=zeros(1,l);     % 预分配内存
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

%% 得到置换加密需要的矩阵h1
%产生Logistic混沌序列p
u0=3.98;     %Logistic参数μ，自定为3.98，可作为密钥
x0=sum(sum(I(1:M/2,:)));     %计算得出Logistic初值x0，可作为密钥
x0=(2*x0)/SUM;
x0=floor(x0*10^4)/10^4;     %保留4位小数
H1=zeros(1,l+1000);        %预分配内存
H1(1)=x0;
for i=1:l+999                 %进行l+999次循环，共得到l+1000点（包括初值）
    H1(i+1)=u0*H1(i)*(1-H1(i));
end
H1=H1(1001:length(H1));            %去除前1000点，获得更好的随机性
[~,h1]=sort(H1,'descend');

%% 得到符号加密需要的矩阵H2
u1=3.892;     %Logistic参数μ，自定为3.892，可作为密钥
x1=sum(sum(I(M/2+1:M,:)));     %计算得出Logistic初值x0，可作为密钥
x1=(2*x1)/SUM;
x1=floor(x1*10^4)/10^4;     %保留4位小数
H2=zeros(1,l+1000);        %预分配内存
H2(1)=x1;
for i=1:l+999                 %进行L+999次循环，共得到L+1000点（包括初值）
    H2(i+1)=u1*H2(i)*(1-H2(i));
end
H2=H2(1001:length(H2));            %去除前1000点，获得更好的随机性
%将H2转换成-1，1序列
H2(H2>0.5)=1;
H2(H2<=0.5)=-1;

%% 重复加密count次
count = 100;     %重复加密count次，可更改，可作为密钥
B1 = A;
B=zeros(1,l);       %预分配内存
for i = 1:count
    for j=1:l
        B(j)=B1(h1(j));     %置换加密
    end
    B1=B.*H2;               %符号加密
end
save ../加密后数据/jiami.mat B1;
%% 不解密 直接解压
% 将一维序列转换为原来每个分块后的二维序列的左上角部分，其余部分补零
GG=zeros(M,N);
for i=1:SUM/t^2
    j=1;
    B2=B1(1+a*(i-1):i*a); % 原来每个分块中取左上角15个元素
    G3=zeros(rate);
    for m=1:rate
        for n=1:rate+1-m
            G3(m,n)=B2(j);
            j=j+1;
        end
    end
    %确定此分块的横纵坐标
    xx=floor(i/(N/t))+1;      %此分块属于第几大行
    yy=mod(i,N/t);            %此分块第几大列
    if yy==0
        xx=xx-1;
        yy=N/t;
    end
    GG((xx-1)*t+1:(xx-1)*t+rate,(yy-1)*t+1:(yy-1)*t+rate)=G3; % 合并所有分块
end

fun = @(block_struct) T'*block_struct.data*T;       %Matlab建议用blockproc函数代替blkproc函数
II = blockproc(GG,[t t],fun);
% II = blkproc(GG,[8,8],'P1*x*P2',T',T);    %重构图像
II_1=im2uint8(II);
imwrite(II_1,'../测试图片和结果/不解密直接解压后的lena.png','png');
figure;imshow(II_1);title('不解密，直接解压后的图片');
figure;imhist(II_1);title('密文图像直方图');

% LL1 = (N*a)/(2*t);
% LL2 = (2*M)/t;
% for i = 1:LL2
%     for j = 1:LL1
%         II(i,j) = B1((i-1)*LL1+j);
%     end
% end
% II_1=II;
% imwrite(II_1,'../测试图片和结果/不解密直接解压后的lena.png','png');
% figure;imshow(II_1);title('不解密，直接解压后的图片');
% figure;imhist(II_1);title('密文图像直方图');
%% 加密后信息熵
%R通道
T2=imhist(II_1);
S2=sum(T2);
xxs2=0;
for i=1:256
    pp2=T2(i)/S2;
    if pp2~=0
        xxs2=xxs2-pp2*log2(pp2);
    end
end

%% 信息输出
disp('加密成功');
disp('密钥：');
disp(['密钥1：μ0=',num2str(u0),'     密钥2：x0=',num2str(x0),'    密钥3：μ1=',num2str(u1),'    密钥4：x1=',num2str(x1)]);
disp(['密钥5：count=',num2str(count),'   密钥6：N=',num2str(N),'   密钥7：M1=',num2str(M1),'   密钥8：N1=',num2str(N1),'   密钥9：rate=',num2str(rate)]);
disp(['原始图片信息熵=',num2str(xxs1),'密文图片信息熵=',num2str(xxs2)]);