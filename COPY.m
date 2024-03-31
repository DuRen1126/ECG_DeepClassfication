clear,clc;
file_names = {'100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '111', '112',...
              '113', '114', '115', '116', '117', '118', '119', '121', '122', '123', '124', '200',...
              '201', '202', '203', '205', '207', '208', '209', '210', '212', '213', '214', '215',...
              '217', '219', '220', '221', '222', '223', '228', '230', '231', '232', '233', '234'};
data = cell(length(file_names), 1); % 创建一个48x1的空元胞数组
for i = 1:length(file_names)
    [signal, Fs, tm] = rdsamp(['database\mitdb\', file_names{i}], []);
    data{i} = signal;
end
plot(tm(1:1000),data{1}(1:1000));title('ECG Signal of Record 100');xlabel('Time(s)');ylabel('Amplitude(mV)');

ANNOT = cell(length(file_names), 1);
ATRTIME = cell(length(file_names), 1);

for j = 1:length(file_names)
    [ann,anntype,subtype]=rdann(['database\mitdb\', file_names{j}],'atr', [],[],[],'N');
    ANNOT{j} = categorical(cellstr(anntype));
    ATRTIME{j} = tm(ann);
end
