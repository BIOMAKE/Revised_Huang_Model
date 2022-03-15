%%% The Corrected Huang's Model %%%
%tbi - time interval between bursts, sec
%tgap - time interval between trains, sec
%Bt - the number of burst per train, positive integer, larger than one
%T - the total number of trains in one stimulation session, positive
%integer
%pc - indicator if it is preceded by a tonic contraction, 1=yes, 0=no
%fpc - indicator if there is a post muscle contraction immediate after TBS
clear
addpath("ExperimentalMeasurements\") %all experimental data stored in this folder


tbi = 0.16; tgap = 8; Bt = 10; T = 20; pc = 1; fpc = 0;

%%%% Check the parameters: ensure no wrong inputs %%%%
if tgap == 0
    Bt = T*Bt;
    T = 1;
end

%%%% Calcium parameters %%%%
if pc == 1
    C = 1; %peak calcium level after a burst
elseif pc == 0
    C = 3;
end
k = 1.2; %decay constant of calcium level after each burst
Rf = 1.55; %proportionality constant of facilitation
fk = 0.005/C; %decay constant of facilitation
Ri = 0.12/C^2; %proportionality constant of inhibition
ik = 0.07/C^2; %decay constant of inhibition
bk = 0.1; %decay constant for the declines between trains in iTBS and imTBS

%%%% STAGE I & II For One Train%%%%
%STAGE I: variables
Cmin = zeros(1,Bt+1);
%STAGE II: variables
Fmin = zeros(1,Bt+1);
Imin = zeros(1,Bt+1);

for n = 1:Bt
    %%% STAGE I: Calcium accumulation (min) %%%
    Cmin(n+1) = C*(1-exp(-n*k*tbi))/(1-exp(-k*tbi))*exp(-k*tbi);
    
    %%% STAGE II: Facilitation/Inhibition substances (min) %%%
    %Facilitation
    Fmin(n+1) = (Rf*(Cmin(n+1)-Cmin(n)) + Fmin(n))*exp(-fk);

    %Inhibition
    Imin(n+1) = (Ri*Cmin(n+1)+Imin(n))*exp(-ik);
end


%%% Train and decay during intertrain interval %%%
% Gap between trains, time step = 0.1s
base_gap = 0.1:0.1:tgap;
% Whole session
session =  1 + T*Bt + length(base_gap)*(T-1); %start point + during burst train + (T-1)during gaps
Cmin_train = zeros(1,session);
Fmin_train = zeros(1,session);
Imin_train = zeros(1,session);
time_axis = zeros(1,session);

% Recursive calculation
cidx = 1;
for t = 1:T
    %Fill the burst session
    Cmin_train(cidx+1:cidx+Bt) = Cmin_train(cidx) + Cmin(2:end);
    Fmin_train(cidx+1:cidx+Bt) = Fmin_train(cidx) + Fmin(2:end);
    Imin_train(cidx+1:cidx+Bt) = Imin_train(cidx) + Imin(2:end);
    for i = cidx+1:cidx+Bt
        time_axis(i) = time_axis(i-1) + 0.04 + tbi;
    end
    %update index
    cidx = cidx + Bt;
    if t ~= T
        %Fill the gap between trains
        Cmin_train(cidx+1:cidx+length(base_gap)) = Cmin_train(cidx)*exp(-base_gap);
        Fmin_train(cidx+1:cidx+length(base_gap)) = Fmin_train(cidx)*exp(-bk*base_gap);
        Imin_train(cidx+1:cidx+length(base_gap)) = Imin_train(cidx)*exp(-bk*base_gap);
        time_axis(cidx+1:cidx+length(base_gap)) = time_axis(cidx) + base_gap;
        %update index
        cidx = cidx + length(base_gap);
    end
end

%%%% STAGE III: After Effect %%%%
%Parameters for Stage III
B = Bt*T; % total number of bursts
h = 3.6;
Inhi_peak = round(1800*B^h/(130^h+B^h)); % peak time for inhibitory after-effect
Faci_peak = round(Inhi_peak/5); % peak time for facilitory after-effec
FinalFaci = Fmin_train(end); %this final value is sampled right after the last train
FinalInhi = Imin_train(end);

if T > 1
    % short version of final value of facilitation
    FC_check = Fmin(end)*(1-exp(-T*bk*tgap))/(1-exp(-bk*tgap));
    % short version of final value of inhibition
    IC_check = Imin(end)*(1-exp(-T*bk*tgap))/(1-exp(-bk*tgap));
else
    FC_check = Fmin(end);
    IC_check = Imin(end);
end

disp(FinalFaci-FC_check)
disp(FinalInhi-IC_check)

%Net summation of facilitation and inhibition
if fpc == 0 
    net = FinalFaci - FinalInhi;
else
    net = FinalFaci;
end

%Parameters for facilitation
h1f=2.5;
h2f=3;
t50oEXC=Faci_peak*0.2;
t50dEXC=Faci_peak*1.1;

%Parameters for inhibition
h1i=3;
h2i=2;
t50oINH=Inhi_peak*0.25;
t50dINH=Inhi_peak*2;

%Aftereffect: this part of code is original Huang's code.
aevalue=[];
time=[];
if or (FinalInhi >= FinalFaci, pc==0)
    for ae=0:4500
        if ae<=Inhi_peak
           I=net*ae^h1i/(t50oINH^h1i+ae^h1i);
           Ip=I;
        else
           I=Ip*t50dINH^h2i/(t50dINH^h2i+(ae-Inhi_peak)^h2i);
        end
    aevalue=[aevalue I];
    time=[time ae];
    end
else
    for ae=0:4500
        if ae<=Faci_peak
            F=net*ae^h1f/(t50oEXC^h1f+ae^h1f);
            Fp=F;
        else
            F=Fp*t50dEXC^h2f/(t50dEXC^h2f+(ae-Faci_peak)^h2f);
        end
    aevalue=[aevalue F];
    time=[time ae];
    end
end



%%
%load file
fileName = {'iTBS600'};
load(fileName{1},'A');

% %Plot figure
% figure('Renderer', 'painters', 'Position', [10 10 400 500])
% %calcium level
% subplot(311)
% plot(time_axis, Cmin_train, 'k');
% ylabel({'Stage I:','Calcium Level',''})
% ylim([0,4])

%Stage II
subplot(2, 1, 1);
hold on
plot(time_axis, Fmin_train, 'b-', 'LineWidth', 2);

plot(time_axis, Imin_train, 'b--', 'LineWidth', 1.5);

legend('Facilitation','Inhibition','Location','southeast','Orientation','horizontal')

ylabel('Stage II','FontSize',15)

ylim([0,12])
hold off;

%Stage III
subplot(2,1,2)
hold on
plot(A.AE(1,:),A.AE(2,:),'Color','k','LineStyle','--','Marker','o')
plot(time,aevalue,'k');
yline(0,'k--')
hold off
legend('Measured','Simulated')
yticks([-10,-5, 0, 5, 10])
ylim([-15,15])
xlim([0,3300])


ylabel('STAGE III:','FontSize',15)
xlabel('Time in Second')


































