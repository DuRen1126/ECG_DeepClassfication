% 清除命令窗口、工作区和图形窗口
clear; clc; close all;

if ~exist('data_denoised', 'var')
    load('data_denoised.mat');
end

% 读取文件名列表
fileID = fopen('database\mitdb\RECORDS.txt','r');% 打开文件
data = textscan(fileID, '%s');% 读取数据
fclose(fileID);% 关闭文件
file_names = data{1};% 将数据存储到数组中

% 初始化标签样本数组
SAMPLE_N = cell(size(file_names));
SAMPLE_L = cell(size(file_names));
SAMPLE_R = cell(size(file_names));
SAMPLE_V = cell(size(file_names));

ATRTIME_N = cell(length(file_names), 1);
ATRTIME_L = cell(length(file_names), 1);
ATRTIME_R = cell(length(file_names), 1);
ATRTIME_V = cell(length(file_names), 1);

for j = 1:length(file_names)
    [ann,anntype] = rdann(['database\mitdb\', file_names{j}],'atr');
    valid_idx_N = anntype == 'N';%valid_idx_N是一个逻辑索引，用于过滤出标签为 N 的注释类型在 anntype 中的位置
    valid_idx_L = anntype == 'L';%valid_idx_L是一个逻辑索引，用于过滤出标签为 L 的注释类型在 anntype 中的位置
    valid_idx_R = anntype == 'R';%valid_idx_R是一个逻辑索引，用于过滤出标签为 R 的注释类型在 anntype 中的位置
    valid_idx_V = anntype == 'V';%valid_idx_V是一个逻辑索引，用于过滤出标签为 V 的注释类型在 anntype 中的位置
    ATRTIME_N{j} = ann(valid_idx_N);
    ATRTIME_L{j} = ann(valid_idx_L);
    ATRTIME_R{j} = ann(valid_idx_R);
    ATRTIME_V{j} = ann(valid_idx_V);
end
% save('ATRTIME_N.mat', 'ATRTIME_N');
% save('ATRTIME_L.mat', 'ATRTIME_L');
% save('ATRTIME_R.mat', 'ATRTIME_R');
% save('ATRTIME_V.mat', 'ATRTIME_V');

% 初始化样本计数器
count_N = 0;
count_L = 0;
count_R = 0;
count_V = 0;

for j = 1:length(file_names)
    for k = 1:length(ATRTIME_N{j})
        idx = ATRTIME_N{j}(k);
        if idx-150 >= 1 && idx+150 <= length(data_denoised{j})
            count_N = count_N + 1;
            SAMPLE_N{count_N} = data_denoised{j}(idx-150:idx+150);
        end
    end
    save('SAMPLE_N.mat', 'SAMPLE_N');

    for k = 1:length(ATRTIME_L{j})
        idx = ATRTIME_L{j}(k);
        if idx-150 >= 1 && idx+150 <= length(data_denoised{j})
            count_L = count_L + 1;
            SAMPLE_L{count_L} = data_denoised{j}(idx-150:idx+150);
        end
    end
    save('SAMPLE_L.mat', 'SAMPLE_L');

    for k = 1:length(ATRTIME_R{j})
        idx = ATRTIME_R{j}(k);
        if idx-150 >= 1 && idx+150 <= length(data_denoised{j})
            count_R = count_R + 1;
            SAMPLE_R{count_R} = data_denoised{j}(idx-150:idx+150);
        end
    end
    save('SAMPLE_R.mat', 'SAMPLE_R');

    for k = 1:length(ATRTIME_V{j})
        idx = ATRTIME_V{j}(k);
        if idx-150 >= 1 && idx+150 <= length(data_denoised{j})
            count_V = count_V + 1;
            SAMPLE_V{count_V} = data_denoised{j}(idx-150:idx+150);
        end
    end
    save('SAMPLE_V.mat', 'SAMPLE_V');
end


figure;
Fs = 360; % 采样频率
t = -150/Fs:1/Fs:150/Fs; % 时间轴
subplot(2,2,1);
plot(t, SAMPLE_N{1},'LineWidth',1);
xlim([-0.4,0.4]); % 横坐标范围
ylim([-1.5, 1.5]); % 纵坐标范围
xlabel('时间 (s)');
ylabel('电压 (mV)');
title('N样本');
subplot(2,2,2);
plot(t, SAMPLE_L{1},'LineWidth',1);
xlim([-0.4,0.4]);
ylim([-1.5, 1.5]);
xlabel('时间 (s)');
ylabel('电压 (mV)');
title('L样本');
subplot(2,2,3);
plot(t, SAMPLE_R{1},'LineWidth',1);
xlim([-0.4,0.4]);
ylim([-1.5, 1.5]);
xlabel('时间 (s)');
ylabel('电压 (mV)');
title('R样本');
subplot(2,2,4);
plot(t, SAMPLE_V{4},'LineWidth',1);
xlim([-0.4,0.4]);
ylim([-1.5, 1.5]);
xlabel('时间 (s)');
ylabel('电压 (mV)');
title('V样本');

subplot(2,2,1)
xlim([-0.246 0.251])
ylim([-0.82 1.04])
subplot(2,2,2)
xlim([-0.41 0.39])
ylim([-1.26 1.74])
subplot(2,2,3)
xlim([-0.41 0.39])
ylim([-2.28 0.63])
subplot(2,2,4)
xlim([-0.41 0.42])
ylim([-1.72 2.48])

% 统计样本数
fprintf('N样本数：%d\n', count_N);
fprintf('L样本数：%d\n', count_L);
fprintf('R样本数：%d\n', count_R);
fprintf('V样本数：%d\n', count_V);

