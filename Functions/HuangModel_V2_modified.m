%%% Huang's Model %%%
function [AF, Fmin_train, Imin_train, time_axis] = HuangModel_V2_modified(X, Bt, T, tbi, tgap, pc, fpc, X_optimal)
%tbi - time interval between bursts, sec
%tgap - time interval between trains, sec
%Bt - the number of burst per train, positive integer, larger than one
%T - the total number of trains in one stimulation session, positive
%integer
%pc - indicator if it is preceded by a tonic contraction, 1=yes, 0=no
if nargin < 5
    pc = 1;
end
%%%% Check the parameters: ensure no wrong inputs %%%%
if tgap == 0
    Bt = T*Bt;
    T = 1;
end


%%%% Calcium parameters %%%%
if pc == 1
    C = X_optimal(1); %peak calcium level after a burst
elseif pc == 0
    C = X_optimal(2);
end
k = X_optimal(11); %decay constant of calcium level after each burst
bk = X_optimal(12); %decay constant for the declines between trains in iTBS and imTBS

Rf = X_optimal(13); %proportionality constant of facilitation
fk = X_optimal(14)/C; %decay constant of facilitation
Ri = X_optimal(15)/C^2; %proportionality constant of inhibition
ik = X_optimal(16)/C^2; %decay constant of inhibition


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
    if n == 1
        Fmin(n+1) = Rf*Cmin(n+1)*exp(-fk);
    else
        Fmin(n+1) = (Rf*(Cmin(n+1)-Cmin(n)) + Fmin(n))*exp(-fk);
    end
    
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


%Net summation of facilitation and inhibition
if fpc == 0
    net = FinalFaci - FinalInhi;
else
    net = FinalFaci;
end

%Parameters for facilitation
h1f=X_optimal(7);
h2f=X_optimal(8);
t50oEXC=Faci_peak*X_optimal(9);
t50dEXC=Faci_peak*X_optimal(10);

%Parameters for inhibition
h1i=X_optimal(3);
h2i=X_optimal(4);
t50oINH=Inhi_peak*X_optimal(5);
t50dINH=Inhi_peak*X_optimal(6);

%after effect
AF = zeros(size(X));
for i = 1:length(X)
    ae = X(i);
    if or (FinalInhi >= FinalFaci, pc==0)
        Ip = net*Inhi_peak^h1i/(t50oINH^h1i+Inhi_peak^h1i);
        if ae <= Inhi_peak
            AF(i) = net*ae^h1i/(t50oINH^h1i+ae^h1i);
        else
            AF(i) = Ip*t50dINH^h2i/(t50dINH^h2i+(ae-Inhi_peak)^h2i);
        end

    else
        Fp = net*Faci_peak^h1f/(t50oEXC^h1f+Faci_peak^h1f);
        if ae <= Faci_peak
            AF(i) = net*ae^h1f/(t50oEXC^h1f+ae^h1f);
        else
            AF(i) = Fp*t50dEXC^h2f/(t50dEXC^h2f+(ae-Faci_peak)^h2f);
        end
    end

end



end





