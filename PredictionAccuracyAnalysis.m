clear
addpath("Functions\")
addpath("ExperimentalMeasurements\")

%load measurement file

%fileName = {'iTBS600', 'imTBS600', 'cTBS600_60min'};

%fileName = {'cTBS300', 'cTBS300_noPC', 'cTBS600_noPC'};

%fileName = {'cTBS300_AC', 'iTBS600_AC'};

fileName = {'cTBS300','cTBS600_60min','imTBS600','iTBS600','cTBS300_noPC','cTBS600_noPC','iTBS600_AC','cTBS300_AC'};

X_optimum = [1, 3, 2.5, 4, 0.2, 1.1, 3, 2, 0.25, 2, 1.2, 0.1, 1.45, 0.005, 0.123, 0.070];



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

new_rmse = sqrt(sum((measure-Modified_predict).^2)/length(measure));
old_rmse = sqrt(sum((measure-Old_predict).^2)/length(measure));

%% correlation coefficient

new_r = corrcoef(measure, Modified_predict);
old_r = corrcoef(measure, Old_predict);

%% r2
R2_new = 1 - sum((measure-Modified_predict).^2)/sum((measure-mean(measure)).^2);
R2_old = 1 - sum((measure-Old_predict).^2)/sum((measure-mean(measure)).^2);

error_new = measure-Modified_predict;
std_new = std(error_new);

error_old = measure-Old_predict;
std_old = std(error_old);



%goodness of fit figure
figure
set(gcf,'unit','centimeters','position',[5,5,15,10]);
set(gcf,'color','w');
set(gcf,'defaultAxesTickLabelInterpreter','latex');

t1 = tiledlayout('flow',"TileSpacing","compact","Padding","compact");

x = 1:3;
y = [R2_new, R2_old; new_rmse,old_rmse; std_new, std_old];
b = bar(x,y);
xname = {'Pseudo R-squared', 'RMSE', 'Error SD'};
xticklabels(xname)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',12)
b(2).FaceColor = '#A2142F';

legend('Revised Model','Initial Model','Location','northwest',...
    'Orientation','vertical','Interpreter','latex','FontSize',12)

exportgraphics(gcf,'Goodness of fit.pdf','ContentType','vector');




%basic setting
%set figure size
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[5,5,20,10]);

t = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

%plot 1
f1 = nexttile;
histogram(error_old,'BinWidth',1,'FaceColor','k')
xlabel(f1,'\textbf{(a)}','Interpreter','latex','FontSize',18)
set(f1,'YLim',[0,35],'XLim',[-10,10])
title(f1,'\textbf{Initial Model}','FontSize',18,'Interpreter','latex')

%plot 2
f2 = nexttile;
h1 = histogram(error_new,'BinWidth',1);
h1.FaceColor = 'k';
set(f2,'YLim',[0,35],'XLim',[-10,10])
title(f2,'\textbf{Revised Model}','FontSize',18,'Interpreter','latex')
xlabel(f2,'\textbf{(b)}','Interpreter','latex','FontSize',18)

ylabel(t,'\textbf{Frequency}','Interpreter','latex','FontSize',18)


%export images
exportgraphics(gcf,'Residual distribution.pdf','ContentType','vector');














