warning off;
% clear,clc,close all;

if ~exist('SAMPLE_L', 'var') || ~exist('SAMPLE_N', 'var') || ~exist('SAMPLE_R', 'var') || ~exist('SAMPLE_V', 'var')
    load('SAMPLE_L.mat');
    load('SAMPLE_N.mat');
    load('SAMPLE_R.mat');
    load('SAMPLE_V.mat');
end

features_L = table; 
features_N = table; 
features_R = table; 
features_V = table; 

for i = 1:length(SAMPLE_L(1:5000))
    
    v =  cell2mat(SAMPLE_L(i));
    %时域特征
    features_L.Mean(i) = mean(v);                         %平均值
    features_L.Std(i) = std(v);                           %标准差
    features_L.Skewness(i) = skewness(v);                 %偏度
    features_L.Kurtosis(i) = kurtosis(v);                 %峭度
    features_L.max(i) = max(v);                           %最大值
    features_L.min(i) = min(v);                           %最小值
    features_L.Peak2Peak(i) = peak2peak(v);               %峰峰值
    features_L.RMS(i) = rms(v);                           %均方根
    features_L.CrestFactor(i) = max(v)/rms(v);            %振幅因数
    features_L.ShapeFactor(i) = rms(v)/mean(abs(v));      %波形因数
    features_L.ImpulseFactor(i) = max(v)/mean(abs(v));    %冲击因数
    features_L.MarginFactor(i) = max(v)/mean(abs(v))^2;   %裕度因数
    features_L.Energy(i) = sum(v.^2);                     %能量
    features_L.Label(i) = 'L';
end

for i = 1:length(SAMPLE_N(1:5000))
    
    v =  cell2mat(SAMPLE_N(i));
    %时域特征
    features_N.Mean(i) = mean(v);                         %平均值
    features_N.Std(i) = std(v);                           %标准差
    features_N.Skewness(i) = skewness(v);                 %偏度
    features_N.Kurtosis(i) = kurtosis(v);                 %峭度
    features_N.max(i) = max(v);                           %最大值
    features_N.min(i) = min(v);                           %最小值
    features_N.Peak2Peak(i) = peak2peak(v);               %峰峰值
    features_N.RMS(i) = rms(v);                           %均方根
    features_N.CrestFactor(i) = max(v)/rms(v);            %振幅因数
    features_N.ShapeFactor(i) = rms(v)/mean(abs(v));      %波形因数
    features_N.ImpulseFactor(i) = max(v)/mean(abs(v));    %冲击因数
    features_N.MarginFactor(i) = max(v)/mean(abs(v))^2;   %裕度因数
    features_N.Energy(i) = sum(v.^2);                     %能量
    features_N.Label(i) = 'N';
end

for i = 1:length(SAMPLE_R(1:5000))
    
    v =  cell2mat(SAMPLE_R(i));
    %时域特征
    features_R.Mean(i) = mean(v);                         %平均值
    features_R.Std(i) = std(v);                           %标准差
    features_R.Skewness(i) = skewness(v);                 %偏度
    features_R.Kurtosis(i) = kurtosis(v);                 %峭度
    features_R.max(i) = max(v);                           %最大值
    features_R.min(i) = min(v);                           %最小值
    features_R.Peak2Peak(i) = peak2peak(v);               %峰峰值
    features_R.RMS(i) = rms(v);                           %均方根
    features_R.CrestFactor(i) = max(v)/rms(v);            %振幅因数
    features_R.ShapeFactor(i) = rms(v)/mean(abs(v));      %波形因数
    features_R.ImpulseFactor(i) = max(v)/mean(abs(v));    %冲击因数
    features_R.MarginFactor(i) = max(v)/mean(abs(v))^2;   %裕度因数
    features_R.Energy(i) = sum(v.^2);                     %能量
    features_R.Label(i) = 'R';
end

for i = 1:length(SAMPLE_V(1:5000))
    
    v =  cell2mat(SAMPLE_V(i));
    %时域特征
    features_V.Mean(i) = mean(v);                         %平均值
    features_V.Std(i) = std(v);                           %标准差
    features_V.Skewness(i) = skewness(v);                 %偏度
    features_V.Kurtosis(i) = kurtosis(v);                 %峭度
    features_V.max(i) = max(v);                           %最大值
    features_V.min(i) = min(v);                           %最小值
    features_V.Peak2Peak(i) = peak2peak(v);               %峰峰值
    features_V.RMS(i) = rms(v);                           %均方根
    features_V.CrestFactor(i) = max(v)/rms(v);            %振幅因数
    features_V.ShapeFactor(i) = rms(v)/mean(abs(v));      %波形因数
    features_V.ImpulseFactor(i) = max(v)/mean(abs(v));    %冲击因数
    features_V.MarginFactor(i) = max(v)/mean(abs(v))^2;   %裕度因数
    features_V.Energy(i) = sum(v.^2);                     %能量
    features_V.Label(i) = 'V';
end

save('features_L.mat', 'features_L');
save('features_N.mat', 'features_N');
save('features_R.mat', 'features_R');
save('features_V.mat', 'features_V');

feature_train = [features_L(1:4500,:); features_N(1:4500,:); features_R(1:4500,:); features_V(1:4500,:)];
feature_test = [features_L(4501:5000,:); features_N(4501:5000,:); features_R(4501:5000,:); features_V(4501:5000,:)];

