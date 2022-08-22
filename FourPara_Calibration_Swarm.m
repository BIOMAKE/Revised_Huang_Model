%Optimisation for the modified Huang's model
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
%{0,1,0,0,1,0,0,1}

%data extraction from dataset
protocols = cell(size(fileName));
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
reg = 0;
sigma = 0;
J = @(x) sumCostFun(protocols, C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk, ...
    x(1), x(2), x(3), x(4)) + reg*(norm(x))^2 + sigma*(norm(L1*x'))^2;



%optimization method
A = [];
b = [];
Aeq = [];
beq = [];


lb = [0.001, 0.001, 0.001, 0.001];
ub = [2, 1, 1, 1];



options = optimoptions('particleswarm', 'PlotFcns','pswplotbestf', 'Display', 'iter');
[X_optimum, fval, exitflag, output] = particleswarm(J,n,lb,ub,options); %optimisation result



%Error criterion
F1 = J(X_optimum)/(116-16+1);



%plot the fitted curve
for i = 1:length(fileName)
    pattern = protocols{i}.pattern;
    pc = protocols{i}.pc;
    [Faci, Inhi] = peakM(pattern, C1, C2, pc, k, bk, ...
        X_optimum(1), X_optimum(2), X_optimum(3), X_optimum(4));
    M = [Faci, Inhi];

    figure
    hold on
    if or (Faci <= Inhi, pc == 0)
        PlotFittingCurve(fileName{i}, protocols{i}, M, [h1in,h2in,k1in,k2in])
    else %facilitation
        PlotFittingCurve(fileName{i}, protocols{i}, M, [h1fa,h2fa,k1fa,k2fa])
    end
    hold off
end

X_optimum = [C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk, X_optimum];

%%
%Cost function
function J = sumCostFun(protocols, C1, C2, h1in, h2in, k1in, k2in, h1fa, h2fa, k1fa, k2fa, k, bk, D1, D2, D3, D4)
%C1 - calcium influx with prior contraction
%C2 - calcium influx without prior contraction
J = 0;
for i = 1:length(protocols)
    pattern = protocols{i}.pattern;
    pc = protocols{i}.pc;
    [Faci, Inhi] = peakM(pattern, C1, C2, pc, k, bk, D1, D2, D3, D4);
    if or (Faci <= Inhi, pc == 0)
        J = J + ModifiedCostFun(protocols{i},C1,C2,h1in,h2in,k1in,k2in,k,bk,D1,D2,D3,D4);
    else
        J = J + ModifiedCostFun(protocols{i},C1,C2,h1fa,h2fa,k1fa,k2fa,k,bk,D1,D2,D3,D4);
    end
end

end























