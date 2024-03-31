if ~exist('ARMA_coeffs_L', 'var') || ~exist('ARMA_coeffs_N', 'var') || ~exist('ARMA_coeffs_R', 'var') || ~exist('ARMA_coeffs_V', 'var')
    load('ARMA_coeffs_L.mat');
    load('ARMA_coeffs_N.mat');
    load('ARMA_coeffs_R.mat');
    load('ARMA_coeffs_V.mat');
end

% 运行下面的代码

% 设置随机种子
rng(1);

% 获取样本总数和特征数
[n_samples, n_features] = size(ARMA_coeffs);

% % 随机划分训练集、验证集和测试集
% train_ratio = 0.7;
% val_ratio = 0.2;
% test_ratio = 0.1;

n_train = round(train_ratio * n_samples);
n_val = round(val_ratio * n_samples);
n_test = n_samples - n_train - n_val;

% 打乱数据集
idx = randperm(n_samples);

% 划分数据集
train_set = ARMA_coeffs(idx(1:n_train), :);
train_labels = labels(idx(1:n_train));

val_set = ARMA_coeffs(idx(n_train+1:n_train+n_val), :);
val_labels = labels(idx(n_train+1:n_train+n_val));

test_set = ARMA_coeffs(idx(n_train+n_val+1:end), :);
test_labels = labels(idx(n_train+n_val+1:end));
