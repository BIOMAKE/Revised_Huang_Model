%plot figures
clear
addpath("Functions\")
addpath("ExperimentalMeasurements\")

%set figure size
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[5,5,35,20]);
%set(gca,'units','centimeters','position',[1,1,19,14]);


%linewidth and marksize
linewidth = 2;
marksize = 8;

legendfont = 15;
labelfont = 20;
titilefont = 20;

%plot data
fileName = {'cTBS300', 'cTBS300_noPC', 'cTBS600_noPC'};

X_optimum = [1, 3, 2.5, 4, 0.2, 1.1, 3, 2, 0.25, 2, 1.2, 0.1, 1.45, 0.005, 0.123, 0.070];

n = length(fileName);

%color map
colorset = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], 'b', 'r', 'k'};

%panel stage II
%p1 = uipanel('Position',[.01 .51 .98 .47]);
t1 = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');


%panel stage III
%p2 = uipanel('Position',[.01 .02 .98 .47]);
%t2 = tiledlayout(p2,1,3,'TileSpacing','Compact','Padding','Compact');




for i = 1:length(fileName)
    %measurement data
    load(fileName{i},'A');
    time = A.AE(1,:)';
    measure = A.AE(2,:)';
    pattern = A.pattern;

    time_new = 0:1:4000;
    %old prediction
    [old_prediction, old_Faci, old_Inhi, old_time] = ...
        HuangModel_Old(time_new,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc);

    %new prediction
    [new_prediction, new_Faci, new_Inhi, new_time] = ...
        HuangModel_V2_modified(time_new,pattern(1),pattern(2),pattern(3),pattern(4),A.pc,A.fpc,X_optimum);


    %plot stage ii
    h(i) = nexttile(i);
    box on
    %ylim([0 14])
    hold on
    %new
    p1 = plot(new_time,new_Faci,'.-','Color',colorset{1}, 'LineWidth',2, 'DisplayName','Facilitation (Revised)');
    p2 = plot(new_time,new_Inhi,':','Color',colorset{1}, 'LineWidth',2, 'DisplayName','Inhibition (Revised)');
    %old
    p3 = plot(old_time,old_Faci, '.-', 'Color', colorset{2}, 'LineWidth',2,'DisplayName','Facilitation (Initial)');
    p4 = plot(old_time,old_Inhi, ':', 'Color', colorset{2}, 'LineWidth',2,'DisplayName','Inhibition (Initial)');



    %plot stage iii
    h(i+n) = nexttile(i+n);
    box on
    hold on
    ylim([-13 13])
    
    colororder({'k','#D95319'})
    
    yyaxis right
    p5 = plot(time,measure,'LineStyle','--','LineWidth',1.5,...
        'Marker','o','MarkerSize', 6, 'DisplayName','Experimental');
    ytickformat('percentage')
    ylim([-13 13])
    yticklabels({'-100%', '-50%', '0%', '50%', '100%'})

    
    yyaxis left
    p6 = plot(time_new,old_prediction, 'LineStyle','--','Color','r',...
        'LineWidth',2,'DisplayName','Simulated (Initial)');
    p7 = plot(time_new,new_prediction, 'LineStyle','-','Color','b',...
        'LineWidth',2,'DisplayName','Simulated (Revised)');
    yline(0,'LineStyle','--','Color','k')

end

%t1 setting
ylabel(h(1),{'\textbf{Stage II}','\textbf{Levels of substances}'},...
    'FontSize',labelfont,'Interpreter','latex')

title(h(1),'\textbf{A. cTBS300}', '\textbf{with prior contraction}',...
    'Interpreter','latex','FontSize',titilefont)
title(h(2),'\textbf{B. cTBS300}', '\textbf{w/o prior contraction}',...
    'Interpreter','latex','FontSize',titilefont)
title(h(3),'\textbf{C. cTBS600}', '\textbf{w/o prior contraction}',...
    'Interpreter','latex','FontSize',titilefont)
lg = legend(h(3),[p1,p2,p3,p4],'Interpreter', 'latex','FontSize',legendfont,...
    'Location','best');
set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.75]));

%t2 setting
ylabel(h(4),{'\textbf{Stage III}','\textbf{After-effect}'},'FontSize',labelfont,...
            'Interpreter','latex','Color','k')

xlabel(t1,'\textbf{Time in Second}','FontSize',labelfont,'Interpreter','latex')

set(h(4),'XLim',[0,1500])
set(h(5),'XLim',[0,1500])

lg = legend(h(6),[p5,p6,p7],'Interpreter', 'latex','FontSize',legendfont);
set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.75]));

%export images
exportgraphics(gcf,'AfterEffect without prior contraction.pdf','ContentType','vector');

















