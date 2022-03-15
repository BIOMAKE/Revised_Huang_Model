%new after-effect model in Huang's Model
function AT_f = AfterEffectFun(protocol,x,M,h1,h2,k1,k2)
%Inputs:
%protocol - structure data of TBS protocol, pattern and experiment records
%x - time
%M - peak after effect
%h1,h2,k1,k2 - parameters to be fitted

%extract data from 
pattern = protocol.pattern;
pc = protocol.pc;
fpc = protocol.fpc;
FinalFaci = M(1);
FinalInhi = M(2);

if fpc == 0
    net = FinalFaci-FinalInhi;
else
    net = FinalFaci;
end


%peak time
Totalpulse = pattern(1)*pattern(2);
h = 3.6;
tpeakInhi = round(1800*Totalpulse^h/(130^h + Totalpulse^h));
tpeakFaci = round(tpeakInhi/5);

%Plot the fitted curve
AT_f = zeros(size(x));
if or (FinalInhi >= FinalFaci, pc == 0)
    Ip = net*tpeakInhi^h1/((k1*tpeakInhi)^h1 + tpeakInhi^h1);
    for i = 1:length(x)
        if x(i) <= tpeakInhi
            AT_f(i) = net*x(i)^h1/((k1*tpeakInhi)^h1 + x(i)^h1);
        else
            AT_f(i) = Ip*(k2*tpeakInhi)^h2/((k2*tpeakInhi)^h2 + (x(i)-tpeakInhi)^h2);
        end
    end
else
    Fp = net*tpeakFaci^h1/((k1*tpeakFaci)^h1 + tpeakFaci^h1);
    for i = 1:length(x)
        if x(i) <= tpeakFaci
            AT_f(i) = net*x(i)^h1/((k1*tpeakFaci)^h1 + x(i)^h1);
        else
            AT_f(i) = Fp*(k2*tpeakFaci)^h2/((k2*tpeakFaci)^h2 + (x(i)-tpeakFaci)^h2);
        end
    end
end

end