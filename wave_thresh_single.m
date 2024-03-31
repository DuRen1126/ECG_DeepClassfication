%% ===============================小波阈值去噪============================= %%

%应用db5作为小波函数进行三层分解
%利用无偏似然估计阈值
%对100.dat的单导联数据进行去噪

clear,clc;
%读取数据文件100.dat并将MLII导联数据存储在变量E1中
[signal, ~, ~] = rdsamp('database\mitdb\100.dat', 1);
M = signal;
E1=M(:,1);   % M为第二篇中对100.dat文件处理后得到的数据矩阵，M(:,1)指MLII导联数据，M(:,2)指V5导联数据
E1=E1';
n1=size(E1);
s1=E1(1:2000);
%小波分解
[C1, L1]=wavedec(E1,3,'db5');
%从C中提取尺度3下的近似小波系数
cA3_1=appcoef(C1,L1,'db5',3);
%从信号C中提取尺度1,2,3下的细节小波系数
cD1_1=detcoef(C1,L1,1);
cD2_1=detcoef(C1,L1,2);
cD3_1=detcoef(C1,L1,3);
%使用stein的无偏似然估计原理进行选择各层的阈值
%cD1_1,cD2_1,cD3_1为各层小波系数
%rigrsure为无偏似然估计的阈值类型
thr1_1=thselect(cD1_1,'rigrsure');
thr2_1=thselect(cD2_1,'rigrsure');
thr3_1=thselect(cD3_1,'rigrsure');
%各层的阈值
TR1=[thr1_1,thr2_1,thr3_1];
%'s'为软阈值，'h'为硬阈值
SORH1='s';

%----------去噪----------
% XC为去噪后信号
% [CXC,LXC]为小波分解结构
% PERF0和PERF2是恢复和压缩的范数百分比
% 'lvd'为允许设置各层的阈值
% 'gbl'为固定阈值
% 3为阈值的长度
[XC1,CXC1,LXC1,PERF0_1,PERF2_1]=wdencmp('lvd',E1,'db5',3,TR1,SORH1);

%----------去噪效果衡量----------
%SNR越大效果越好，MSE越小越好
%选取信号的长度
N1=n1(2);
x1=E1;
y1=XC1;
F1=0;
MM1=0;
for ii=1:N1
    m1(ii)=(x1(ii)-y1(ii))^2;
    t1(ii)=y1(ii)^2;
    f1(ii)=t1(ii)/m1(ii);
    F1=F1+f1(ii);
    MM1=MM1+m1(ii);
end
SNR1=10*log10(F1);
MSE1=MM1/N1;
SM1=SNR1/MSE1;

%打印结果
fprintf('SNR1=%f\n', SNR1);
fprintf('MSE1=%f\n', MSE1);

