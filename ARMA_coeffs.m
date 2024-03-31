clear,clc;

%[EstMdl,EstParamCov,logL,info] = estimate(mdl,y)

% 其中，输入参数 mdl 表示ARIMA模型对象，y 表示待拟合的时间序列数据，输出参数包括：
% 
% EstMdl：表示经过拟合得到的ARIMA模型对象；
% EstParamCov：表示模型参数的协方差矩阵；
% logL：表示模型拟合后的对数似然函数值；
% info：是一个结构体，包含拟合过程中的一些信息，如算法使用的优化器、是否收敛等。
% 其中，EstMdl 和 EstParamCov 是后续分析ARIMA模型的关键输出参数。
% EstMdl 包含了拟合得到的ARIMA模型的所有参数信息，包括阶数、系数、方差等；
% 而 EstParamCov 是一个协方差矩阵，用于衡量拟合得到的ARIMA模型参数的精度。

% 判断工作区是否有心搏样本数据，有则跳过，无则加载
if ~exist('SAMPLE_L', 'var') || ~exist('SAMPLE_N', 'var') || ~exist('SAMPLE_R', 'var') || ~exist('SAMPLE_V', 'var')
    load('SAMPLE_L.mat');
    load('SAMPLE_N.mat');
    load('SAMPLE_R.mat');
    load('SAMPLE_V.mat');
end

%判断工作区是否有心搏样本的ARMA参数模型，有则跳过，无则加载
% if ~exist('ARMA_coeffs_L', 'var') || ~exist('ARMA_coeffs_N', 'var') || ~exist('ARMA_coeffs_R', 'var') || ~exist('ARMA_coeffs_V', 'var')
%     load('ARMA_coeffs_L.mat');
%     load('ARMA_coeffs_N.mat');
%     load('ARMA_coeffs_R.mat');
%     load('ARMA_coeffs_V.mat');
% end

%判断工作区是否有心搏样本的ARMA参数模型，有则跳过，无则加载
if ~exist('ARMA_coeffs_L', 'var') || ~exist('ARMA_coeffs_N', 'var') || ~exist('ARMA_coeffs_R', 'var') || ~exist('ARMA_coeffs_V', 'var')

% 初始化ARMA模型参数
p = 4; % AR阶数
q = 2; % MA阶数

% 初始化ARMA模型系数
ARMA_coeffs_L = zeros(5000, p+q+1);
ARMA_coeffs_N = zeros(5000, p+q+1);
ARMA_coeffs_R = zeros(5000, p+q+1);
ARMA_coeffs_V = zeros(5000, p+q+1);

% 对每个L型心搏样本进行ARMA模型处理
for i = 1:length(SAMPLE_L(1:5000))
    % 建立ARMA模型
    mdl = arima(p, 0, q);
    
    % 估计ARMA模型
    [est_mdl, ~, ~, ~] = estimate(mdl, SAMPLE_L{i});

    % 保存ARMA模型系数
    ARMA_coeffs_L(i, :) = [est_mdl.Constant, est_mdl.AR{1}, est_mdl.AR{2}, est_mdl.AR{3}, est_mdl.AR{4}, est_mdl.MA{1}, est_mdl.MA{2}];
end
save('ARMA_coeffs_L.mat', 'ARMA_coeffs_L');

% 对每个N型心搏样本进行ARMA模型处理
for i = 1:length(SAMPLE_N(1:5000))
    % 建立ARMA模型
    mdl = arima(p, 0, q);
    
    % 估计ARMA模型
    [est_mdl, ~, ~, ~] = estimate(mdl, SAMPLE_N{i});

    % 保存ARMA模型系数
    ARMA_coeffs_N(i, :) = [est_mdl.Constant, est_mdl.AR{1}, est_mdl.AR{2}, est_mdl.AR{3}, est_mdl.AR{4}, est_mdl.MA{1}, est_mdl.MA{2}];
end
save('ARMA_coeffs_N.mat', 'ARMA_coeffs_N');

% 对每个R型心搏样本进行ARMA模型处理
for i = 1:length(SAMPLE_R(1:5000))
    % 建立ARMA模型
    mdl = arima(p, 0, q);
    
    % 估计ARMA模型
    [est_mdl, ~, ~, ~] = estimate(mdl, SAMPLE_R{i});

    % 保存ARMA模型系数
    ARMA_coeffs_R(i, :) = [est_mdl.Constant, est_mdl.AR{1}, est_mdl.AR{2}, est_mdl.AR{3}, est_mdl.AR{4}, est_mdl.MA{1}, est_mdl.MA{2}];
end
save('ARMA_coeffs_R.mat', 'ARMA_coeffs_R');

% 对每个V型心搏样本进行ARMA模型处理
for i = 1:length(SAMPLE_V(1:5000))
    % 建立ARMA模型
    mdl = arima(p, 0, q);
    
    % 估计ARMA模型
    [est_mdl, ~, ~, ~] = estimate(mdl, SAMPLE_V{i});

    % 保存ARMA模型系数
    ARMA_coeffs_V(i, :) = [est_mdl.Constant, est_mdl.AR{1}, est_mdl.AR{2}, est_mdl.AR{3}, est_mdl.AR{4}, est_mdl.MA{1}, est_mdl.MA{2}];
end
save('ARMA_coeffs_V.mat', 'ARMA_coeffs_V');

end
