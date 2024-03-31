if ~exist('data_original', 'var')
    load('data_original.mat');
end
Fs = 360;
% 初始化参数
wavelet = 'db4'; % 选择小波基
SNR = 5; % 选择信噪比阈值
L = length(data_original{1}); % 信号长度
N = 2^nextpow2(L); % FFT点数
f = linspace(0, Fs/2, N/2+1); % 频率坐标
T = L/Fs; % 信号持续时间
t = linspace(0, T, L); % 时间坐标

% 对每条信号进行小波去噪
for i = 1:length(data_original)
    % 分解信号
    [C, L] = wavedec(data_original{i}, wmaxlev(length(data_original{i}), wavelet), wavelet);
    % 计算噪声方差
    noise_std = median(abs(C))/0.6745;
    noise_var = noise_std^2;
    % 计算信噪比
    P_signal = sum(C(1:L(1)).^2);
    SNR_i = 10*log10(P_signal/noise_var);
    % 确定分解层数
    level = wmaxlev(length(data_original{i}), wavelet);
    while SNR_i < SNR && level > 1
        level = level - 1;
        C(L(level)+1:end) = 0;
        SNR_i = 10*log10(sum(C(1:L(level)).^2)/noise_var);
    end
    % 基于软阈值进行去噪
    for j = 1:level
        ind = L(j)+1:L(j+1);
        T = noise_std*sqrt(2*log(length(ind)));
        C(ind) = sign(C(ind)).*max(abs(C(ind))-T, 0);
    end
    % 重构信号
    data_denoised{i} = waverec(C, L, wavelet);

    %----------去噪效果衡量----------
%SNR越大效果越好，MSE越小越好
%选取信号的长度
N1=650000;
x1=data_original{i};
y1=data_denoised{i};
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
fprintf('Signal %d: SNR = %f, MSE = %f\n', i, snr_val, mse_val);
end


% 将去噪后的心电信号保存为mat文件
% save('data_denoised.mat', 'data_denoised');

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
