%% 基于DCT变换数字图像压缩解密系统
%   @author:沈洋
%   @date:2018.03.22
%-------------------------------------------------------------------------------------------------------%
load ../加密后数据/jiami.mat;        %导入加密序列，存储时变量名为B1
AA=B1;
% 重复解密count次
count=99;       %密钥，重复加密次数，与加密时保持相同
t=8;            %分块大小
rate=5;     %压缩矩阵1存在的行列数，密钥，与加密时保持相同
a=sum(1:rate);      %每一个子块保留了a个数据
l=length(AA);       %获得加密序列长度
N=512;          %密钥，与加密时保持相同
M=((l/a)*t^2)/N;    %根据N推算得到M
SUM=M*N;
%% 得到置换加密需要的矩阵h1
%产生Logistic混沌序列p
u0=3.98;     %Logistic参数μ，密钥，与加密时保持相同
x0=0.3951;     %Logistic初值x0，密钥，与加密时保持相同
H1=zeros(1,l+1000);        %预分配内存
H1(1)=x0;
for i=1:l+999                 %进行l+999次循环，共得到l+1000点（包括初值）
    H1(i+1)=u0*H1(i)*(1-H1(i));
end
H1=H1(1001:length(H1));            %去除前1000点，获得更好的随机性
[~,h1]=sort(H1,'descend');

%% 得到符号加密需要的矩阵H2
u1=3.892;       %Logistic参数μ，密钥，与加密是保持相同
x1=0.347;      %Logistic初值x0，密钥，与加密时保持相同
H2=zeros(1,l+1000);        %预分配内存
H2(1)=x1;
for i=1:l+999                 %进行L+999次循环，共得到L+1000点（包括初值）
    H2(i+1)=u1*H2(i)*(1-H2(i));
end
H2=H2(1001:length(H2));            %去除前1000点，获得更好的随机性
%将H2转换成-1，1序列
H2(H2>0.5)=1;
H2(H2<0.5)=-1;

%% 开始解密
for i=1:count
    B3=AA.*H2;          %符号解密
    for j=1:l
        AA(h1(j))=B3(j);        %置换解密
    end
end

%将一维序列转化成二维序列，其余部分补零
GG=zeros(M,N);
for i=1:SUM/t^2
    j=1;
    B4=AA(1+a*(i-1):i*a);     % 原来每个分块中取左上角15个元素
    %将15个元素重新拼到分块的左上角
    G3=zeros(rate);
    for m=1:rate
        for n=1:rate+1-m
            G3(m,n)=B4(j);
            j=j+1;
        end
    end
    %确定此分块的横纵坐标
    xx=floor(i/(N/t))+1;      %此分块属于第几大行
    yy=mod(i,(N/t));          %此分块属于第几大列
    if yy==0
        xx=xx-1;
        yy=(N/t);
    end
    GG((xx-1)*t+1:(xx-1)*t+rate,(yy-1)*t+1:(yy-1)*t+rate)=G3;     %将每个分块合并在一起，合成完整的二维图像，其余部分为0
end
T=dctmtx(8);        %8阶DCT变换矩阵
fun = @(block_struct) T'*block_struct.data*T;       %Matlab建议用blockproc函数代替blkproc函数
I_jiemi = blockproc(GG,[t t],fun);              %DCT逆变换，重构图像
% III  = blkproc(GG,[8,8],'P1*x*P2',T',T);    %DCT逆变换，重构图像

%% 去除加密时补上的零
M1=0;   %密钥，与加密时保持相同
N1=0;   %密钥，与加密时保持相同
if M1~=2
    I_jiemi=I_jiemi(1:M-t+M1,:,:);
end
if N1~=0
    I_jiemi=I_jiemi(:,1:N-t+N1,:);
end
imwrite(I_jiemi,'../测试图片和结果/解密后的lena.png','png');
disp('您输入的解密密钥为：');
disp(['密钥1：μ0=',num2str(u0),'     密钥2：x0=',num2str(x0),'    密钥3：μ1=',num2str(u1),'    密钥4：x1=',num2str(x1)]);
disp(['密钥5：count=',num2str(count),'   密钥6：N=',num2str(N),'   密钥7：M1=',num2str(M1),'   密钥8：N1=',num2str(N1),'   密钥9：rate=',num2str(rate)]);
disp('解密完成'); 
figure;imshow(I_jiemi);title(['解密解压图片，压缩率： ',num2str(64),':',num2str(a)]);