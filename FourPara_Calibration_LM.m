%Optimisation for the modified Huang's model by using Levenberg-Marquardt
%Algorithm

%parameters: C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk
%C1 - calcium influx after a burst with prior contraction
%C2 - calcium influx after a burst without prior contraction
%h1in, h2in, k1in, k2in - parameters of inhibitory curve M(t)
%h1fa, h2fa, k1fa, k2fa - parameters of facilitatory curve M(t)
%k - decay time constant of calcium level
%bk - decay time constant of inhibitory/facilitatory substances
%D1 - proportionality constant of facilitatory susbatnce, Rf = D1
%D2 - decay time constant of facilitatory substance, fk = D2/C
%D3 - proportionality constant of inhibitory substance, Ri = D3/C^2
%D4 - decay time constant of inhibitory substance, ik = D4/C^2

clear
addpath("ExperimentalMeasurements\") %all experimental data stored in this folder
addpath("Functions\") %all sub-functions used in this main code

%list of sets of measurements
fileName = {'imTBS600','iTBS600','cTBS300','cTBS600_60min','cTBS300_noPC','cTBS600_noPC','cTBS300_AC','iTBS600_AC'};

%data extraction from dataset
protocols = cell(size(fileName));
xdata = [];
ydata = [];
for i = 1:length(fileName)
    load(fileName{i},'A');
    protocols{i} = A;
end


%initial point for optimisation
C1 = 1; C2 = 3; 

h1in = 2.5; h2in = 4; k1in = 0.2; k2in = 1.1;

h1fa = 3; h2fa = 2; k1fa = 0.25; k2fa = 2;

k = 1.2; bk = 0.1;

D1 = 1; D2 = 0.01; D3 = 0.14; D4 = 0.055; 


%Tikhonov matrix - smooth regularization - smoothness of parameters
n = 4;
L = diag(-ones(1,n-2),2) + diag(2*ones(1,n-1),1) + diag(-ones(1,n));
L = L(1:n-2,:);

%first derivative
L1 = diag(ones(1,n)) + diag(-ones(1,n-1),1);
L1 = L1(1:n-1,:);



%define the cost function
sigma = 0;  %smooth coefficient
mu = 1; %
J = @(x) [sumCostFun(protocols,C1,C2,h1in,h2in,k1in,k2in,h1fa,h2fa,k1fa,k2fa,k,bk, ...
    x(1), x(2), x(3), x(4)); sqrt(sigma)*norm(L*x'); sqrt(mu)*norm(x)];



%optimization method

lb = [];
ub = [];

%lb = [0.001, 0.001, 0.001, 0.001];
%ub = [10, 10, 10, 10];

x0 = [D1, D2, D3, D4]; %start point


options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter-detailed','PlotFcns', 'optimplotresnorm');


[X_optimum,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(J,x0,lb,ub,options);  %optimisation result


%Error criterion
F1 = J(X_optimum)'*J(X_optimum)/(116-16+1);

%confidence interval
ci = nlparci(X_optimum,residual,'jacobian',jacobian);



%plot the fitted curve
% for i = 1:length(fileName)
%     pattern = protocols{i}.pattern;
%     pc = protocols{i}.pc;
%     [Faci, Inhi] = peakM(pattern, C1, C2, pc, k, bk, ...
%         X_optimum(1), X_optimum(2), X_optimum(3), X_optimum(4));
%     M = [Faci, Inhi];
% 
%     figure
%     hold on
%     if or (Faci <= Inhi, pc == 0)
%         PlotFittingCurve(fileName{i}, protocols{i}, M, [h1in,h2in,k1in,k2in])
%     else %facilitation
%         PlotFittingCurve(fileName{i}, protocols{i}, M, [h1fa,h2fa,k1fa,k2fa])
%     end
%     hold off
% end

X_optimum = [C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk, X_optimum];

%%
%Cost function
function J = sumCostFun(protocols, C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk, D1, D2, D3, D4)
%C1 - calcium influx with prior contraction
%C2 - calcium influx without prior contraction
    J = [];
    for i = 1:length(protocols)
        J_mid = [];
        pattern = protocols{i}.pattern;
        pc = protocols{i}.pc;
        tdata = protocols{i}.AE(1,:);
        ydata = protocols{i}.AE(2,:);
    
        [Faci, Inhi] = peakM(pattern, C1, C2, pc, k, bk, D1, D2, D3, D4);
    
        Mnet = [Faci, Inhi];
    
        for j = 1:length(tdata)
            if or (Faci <= Inhi, pc == 0)
                f = AfterEffectFun(protocols{i}, tdata(j), Mnet, h1in, h2in, k1in, k2in) -  ydata(j);
            else
                f = AfterEffectFun(protocols{i}, tdata(j), Mnet, h1fa, h2fa, k1fa, k2fa) -  ydata(j); 
            end
            J_mid = [J_mid; f];
        end
    
        J = [J; J_mid];
    end
    
    
end

























