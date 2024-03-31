% clear,clc;

if ~exist('ARMA_coeffs_L', 'var') || ~exist('ARMA_coeffs_N', 'var') || ~exist('ARMA_coeffs_R', 'var') || ~exist('ARMA_coeffs_V', 'var')
    load('ARMA_coeffs_L.mat');
    load('ARMA_coeffs_N.mat');
    load('ARMA_coeffs_R.mat');
    load('ARMA_coeffs_V.mat');
end

% 绘制混淆矩阵
confusion_val = confusionmat(val_labels, val_pred);
confusion_test = confusionmat(test_labels, test_pred);

% 绘制热力图
figure();
subplot(1,1,1);
h_val = heatmap({'L', 'N', 'R', 'V'}, {'L', 'N', 'R', 'V'}, confusion_val);
xlabel('Predicted labels');
ylabel('True labels');
title('Confusion Matrix (Validation Set)');

figure();
subplot(1,1,1);
h_text = heatmap({'L', 'N', 'R', 'V'}, {'L', 'N', 'R', 'V'}, confusion_test);
xlabel('Predicted labels');
ylabel('True labels');
title('Confusion Matrix (Test Set)');

%初始化acc矩阵
acc_val = zeros(size(confusion_val));
acc_test = zeros(size(confusion_test));

for i = 1:length(confusion_val(1,:))%列循环
    for j = 1:length(confusion_val(:,1))%行循环
        acc_val(i, j) = confusion_val(i, j) / sum(confusion_val(:,i));
    end
end

figure();
subplot(1,1,1);
val_acc = heatmap({'L', 'N', 'R', 'V'}, {'L', 'N', 'R', 'V'}, acc_val);
xlabel('Predicted labels');
ylabel('True labels');
title('Confusion Matrix (Validation Set)');


for i = 1:length(confusion_test(1,:))%列循环
    for j = 1:length(confusion_test(:,1))%行循环
        acc_val(i, j) = confusion_test(i, j) / sum(confusion_test(:,i));
    end
end
figure();
subplot(1,1,1);
test_acc = heatmap({'L', 'N', 'R', 'V'}, {'L', 'N', 'R', 'V'}, acc_val);
xlabel('Predicted labels');
ylabel('True labels');
title('Confusion Matrix (Test Set)');
