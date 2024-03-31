%%%%%%%%%%%%%%%%%%心电信号降噪
%%%%%%%%%%%%%%%Birge-Massart策略阈值降噪
%基于小波变换的心电信号的降噪
% clear, clc, close all;

addpath 'C:\Users\Du\Documents\MATLAB\Examples\R2022b\data_download\mcode'

%判断工作区是否有data_original.mat
if ~exist('data_original', 'var')
    load('data_original.mat');
end

for i = 1:37
wavename='db5';
level=4;
[c,l]=wavedec(data_original{i},level,wavename);
alpha=1.5;
sorh='h';
[thr,nkeep]=wdcbm(c,l,alpha);%使用硬阈值给信号降噪
[xc,cxc,lxc,perf0,perfl2]=wdencmp('lvd',c,l,wavename,level,thr,sorh);
data_denoised{i}=xc;
end

save('data_denoised.mat', 'data_denoised');

figure(1);
subplot(2,1,1);
plot(tm(1:1000),data_original{1}(1:1000));
title('原始心电信号');xlabel('时间(s)');ylabel('电压(mV');

subplot(2,1,2);
plot(tm(1:1000),data_denoised{1}(1:1000));
title('Birge-Massart策略阈值降噪后的ECG信号（wname=db5 level=4）');
xlabel('时间(s)');ylabel('电压(mV');
