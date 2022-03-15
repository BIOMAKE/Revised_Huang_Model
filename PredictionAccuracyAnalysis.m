clear
addpath("Functions\")
addpath("ExperimentalMeasurements\")

%load measurement file

%fileName = {'iTBS600', 'imTBS600', 'cTBS600_60min'};

%fileName = {'cTBS300', 'cTBS300_noPC', 'cTBS600_noPC'};

fileName = {'cTBS300_AC', 'iTBS600_AC'};

%fileName = {'cTBS300','cTBS600_60min','imTBS600','iTBS600','cTBS300_noPC','cTBS600_noPC','iTBS600_AC','cTBS300_AC'};

X_optimum = [1, 3, 2.5, 4, 0.2, 1.1, 3, 2, 0.25, 2, 1.2, 0.1, 1.44, 0.005, 0.127, 0.07];



measure = [];
Modified_predict = [];
Old_predict = [];
for i = 1:length(fileName)
    load(fileName{i},'A');
    time = A.AE(1,:);
    measure = [measure, A.AE(2,:)];
    pattern = A.pattern;
    %prediction
    Modified_predict = [Modified_predict, HuangModel_V2_modified(time,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc,X_optimum)];
    Old_predict = [Old_predict, HuangModel_Old(time,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc)];
end
%% root mean squared error

new_rmse = sqrt(sum((measure-Modified_predict).^2)/length(measure))
old_rmse = sqrt(sum((measure-Old_predict).^2)/length(measure))

%% correlation coefficient

new_r = corrcoef(measure, Modified_predict)
old_r = corrcoef(measure, Old_predict)


%% r2
R2_new = 1 - sum((measure-Modified_predict).^2)/sum((measure-mean(measure)).^2)
R2_old = 1 - sum((measure-Old_predict).^2)/sum((measure-mean(measure)).^2)



error_new = measure-Modified_predict;
std(error_new)
figure
histogram(error_new)
grid on
xlabel('Prediction Error','FontSize',15,'FontWeight','bold')
ylabel('Frequency',FontSize=15,FontWeight='bold')
title('Revised Model',FontSize=15)


error_old = measure-Old_predict;
std(error_old)
figure
histogram(error_old)
grid on
xlabel('Prediction Error','FontSize',15,'FontWeight','bold')
ylabel('Frequency',FontSize=15,FontWeight='bold')
title('Initial Model',FontSize=15)


%%






