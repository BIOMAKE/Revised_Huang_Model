clear
addpath("Functions\")
addpath("ExperimentalMeasurements\")

%load measurement file

fileName = {'iTBS600', 'imTBS600', 'cTBS600_60min'};

%fileName = {'cTBS300', 'cTBS300_noPC', 'cTBS600_noPC'};

%fileName = {'cTBS300_AC', 'iTBS600_AC'};

%fileName = {'cTBS300','cTBS600_60min','imTBS600','iTBS600'};
%fileName = {'cTBS300_noPC','cTBS600_noPC'};
%fileName = {'iTBS600_AC','cTBS300_AC'};

%fileName = {'cTBS300','cTBS600_60min','imTBS600','iTBS600','cTBS300_noPC','cTBS600_noPC','iTBS600_AC','cTBS300_AC'};

X_optimum = [1, 3, 2.5, 4, 0.2, 1.1, 3, 2, 0.25, 2, 1.2, 0.1, 1.44, 0.005, 0.127, 0.07];

n = length(fileName);

figure
for i = 1:length(fileName)
    %measurement data
    load(fileName{i},'A');
    time = A.AE(1,:)';
    measure = A.AE(2,:)';
    pattern = A.pattern;

    time_new = 0:1:4000;
    %old prediction
    [old_prediction, old_Faci, old_Inhi, old_time] = HuangModel_Old(time_new,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc);

    %new prediction
    [new_prediction, new_Faci, new_Inhi, new_time] = HuangModel_V2_modified(time_new,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc,X_optimum);



    %plot stage ii
    subplot(2,n,i)
    hold on
    plot(new_time,new_Faci,'.-','Color',[0.6350 0.0780 0.1840], 'LineWidth',1)
    plot(new_time,new_Inhi,':','Color',[0.6350 0.0780 0.1840], 'LineWidth',2)

    plot(old_time,old_Faci, '.-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth',1)
    plot(old_time,old_Inhi, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth',2)

    if i == length(fileName)
        lgd = legend('Facilitation (Revised)','Inhibition (Revised)', 'Facilitation (Initial)', 'Inhibition (Initial)',...
            'Orientation','horizontal','Location','bestoutside');
        lgd.FontSize = 12;
        hold off
    end



    %plot stage iii
    subplot(2,n,i+n)
    hold on
    ylim([-13 13])
    title(fileName{i},'FontSize',12)
    yyaxis right
    p1 = plot(time,measure,'LineStyle','--','Color','k','LineWidth',1,'Marker','o','DisplayName','Experimental');

    ytickformat('percentage')
    ylim([-13 13])
    yticklabels({'-100%', '-50%', '0%', '50%', '100%'})
    yyaxis left
    p2 = plot(time_new,old_prediction, 'LineStyle','--','Color','r','LineWidth',2,'DisplayName','Prediction (Initial)');
    p3 = plot(time_new,new_prediction, 'LineStyle','-','Color','b','LineWidth',2,'DisplayName','Prediction (Revised)');
    yline(0,'LineStyle','--','Color','k')
    %set(gca,'YTickLabel',{'-100%', '-50%', '0%', '50%', '100%'} )
    
    if i == 1
        ylabel('Normalized after effect','FontSize',12)
    end

    if i == length(fileName)
        legend([p1,p2,p3])
        hold off
    end





end






