%For fitting the after-effects curve of Huang's Model
function [FinalFaci, FinalInhi] = peakM(pattern,C1,C2,pc,k,bk,Rf,fk,Ri,ik)
%Bt - the number of bursts in one train
%T - the number of total trains
%tbi - interval between bursts
%tgap - interval between trains
%C - peak calcium concentration after a burst
%%
Bt = pattern(1); T = pattern(2);
tbi = pattern(3); tgap = pattern(4);

if pc == 1
    C = C1; %prior contraction
else
    C = C2; %without prior contraction
end

%
fk = fk/(C); %decay constant of facilitation
Ri = Ri/(C^2); %proportionality constant of inhibition
ik = ik/(C^2); %decay constant of inhibition


%%
%%%% STAGE I & II For One Train%%%%
%STAGE I: variables
Cmin = zeros(1,Bt+1);

%STAGE II: variables
Fmin = 0; % Min facilitation substance
Imin = 0; % Min inhibition substance

for n = 1:Bt
    %%% STAGE I: Calcium accumulation (min) %%%
    Cmin(n+1) = C*(1-exp(-n*k*tbi))/(1-exp(-k*tbi))*exp(-k*tbi);
    
    %%% STAGE II: Facilitation/Inhibition substances (min) %%%
    %Facilitation
    if n == 1
        Fmin = Rf*Cmin(n+1)*exp(-fk);
    else
        Fmin = (Rf*(Cmin(n+1)-Cmin(n)) + Fmin)*exp(-fk);
    end
   
    %Inhibition
    Imin = (Ri*Cmin(n+1)+Imin)*exp(-ik);
end


%%
%%%% STAGE III: After Effect %%%%
%Parameters for Stage III
if T > 1
    % short version of final value of facilitation
    FinalFaci = Fmin(end)*(1-exp(-T*bk*tgap))/(1-exp(-bk*tgap));
    % short version of final value of inhibition
    FinalInhi = Imin(end)*(1-exp(-T*bk*tgap))/(1-exp(-bk*tgap));
else
    FinalFaci = Fmin(end);
    FinalInhi = Imin(end);
end



end








