%The cost function for one set of measurements
function J = ModifiedCostFun(protocol,C1,C2,h1,h2,k1,k2,k,bk,Rf,fk,Ri,ik)
%Inputs:
%C1 - calcium influx with prior contraction
%C2 - calcium influx without prior contraction
%protocol - TBS protocol
%h1,h2,k1,k2 - parameters to be fitted

%extract data - after effect record
pattern = protocol.pattern;
AE = protocol.AE;
pc = protocol.pc; %pc - with/without prior contraction
x = AE(1,:); %time
y = AE(2,:); %normalized MEP

%calculation of the net effect Mnet
M = zeros(1,2);
[M(1),M(2)] = peakM(pattern,C1,C2,pc,k,bk,Rf,fk,Ri,ik);


%calculation of the corresponding aftereffects
AE_f = AfterEffectFun(protocol,x,M,h1,h2,k1,k2);

%cost for a single set of measurements
J = norm(AE_f-y)^2;

end