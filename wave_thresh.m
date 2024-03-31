%% ===============================小波阈值去噪============================= %%

clear,clc,close all;

fileID = fopen('database\mitdb\RECORDS.txt','r');% 打开文件
data = textscan(fileID, '%s');% 读取数据
fclose(fileID);% 关闭文件
file_names = data{1};% 将数据存储到数组中
% data_original = cell(length(file_names), 1);% 创建一个37x1的空元胞数组

% % 读取数据
% for i = 1:length(file_names)
%     [signal, Fs, tm] = rdsamp(['database\mitdb\', file_names{i}], 1);
%     data_original{i} = signal;
% end

%判断工作区是否有data_original.mat
if ~exist('data_original', 'var')
    load('data_original.mat');
end

%应用db5作为小波函数进行三层分解
%利用无偏似然估计阈值
%对100.dat的单导联数据进行去噪
E1 = data_original;
E1 = E1';
data_denoised = cell(length(file_names), 1);
for j = 1:length(file_names)
n1=size(E1{j});
s1=E1{j}(1:2000);
%小波分解
[C1, L1]=wavedec(E1{j},3,'db5');
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
[XC1,CXC1,LXC1,PERF0_1,PERF2_1]=wdencmp('lvd',E1{j},'db5',3,TR1,SORH1);
data_denoised{j} = XC1;


%----------去噪效果衡量----------
%SNR越大效果越好，MSE越小越好
%选取信号的长度
N1=n1(1);
x1=E1{j};
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
snr_val=10*log10(F1);
mse_val=MM1/N1;
sm_val=snr_val/mse_val;

% 输出结果
%SNR越大越好，MSE越小越好，SM越大越好
fprintf('Signal %d: SNR = %f, MSE = %f, SM = %f\n', j, snr_val, mse_val, sm_val);
end

% 将去噪后的心电信号保存为mat文件
save('data_denoised.mat', 'data_denoised');

% 随机选取三组信号进行比较
figure;
for i = 1:3
    index = randi(length(data_original), 1); % 随机选取一条信号

    % 截取一部分特征，至少包含两个完整的心动周期
    start = 2000; % 起始点
    N = 2000; % 采样点数
    seg_original = data_original{index}(start:start+N-1);
    seg_denoised = data_denoised{index}(start:start+N-1);
     
    % 绘制原始信号和去噪后的信号
    subplot(3,2,2*i-1);
    plot(linspace(0, (N-1)/360, N), seg_original);
    title(sprintf('Original signal %d', index));
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    xlim([0, N/360]);
    ylim([-1.5, 1.5]);
    grid on;
    
    subplot(3,2,2*i);
    plot(linspace(0, (N-1)/360, N), seg_denoised);
    title(sprintf('Denoised signal %d', index));
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    xlim([0, N/360]);
    ylim([-1.5, 1.5]);
    grid on;
end
