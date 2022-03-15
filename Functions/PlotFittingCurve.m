%plot fiiting curve and the original one
function PlotFittingCurve(name, protocol, M, X_optimum)

AE = protocol.AE;
plot(AE(1,:),AE(2,:),'o--','DisplayName',[name, ' Record'])


t = 0:0.1:4500;
plot(t,AfterEffectFun(protocol,t,M,X_optimum(1),X_optimum(2),X_optimum(3),X_optimum(4)),'--','DisplayName',name)


xlabel('t in second')
ylabel('Normalized MEP')
yticks([-10 -5 0 5 10])
axis([0 inf -15 15])
legend show

end