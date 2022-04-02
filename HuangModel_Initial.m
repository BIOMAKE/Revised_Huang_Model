%%%%%               TBS Theoretical Model              %%%%%

clear

burst=25; %burst number
train=8;   %train number
gap=10;     %gap between trains (sec)
pc=1;      %previous tonic contraction 1=yes; 0=no
fbc=0;     %followed by contraction 1=yes; 0=no
bk=.1;     %decay constant of facilitation during gaps


%%%%% Calcium parameters %%%%%
CAmaxmin=[0];
CAx=[0];
C=1;
k=1.2;
t=0.16;
kt=-k*t;

%%%%% facilitation/inhibition substance parameters %%%%%
fmin=[0];
fmaxmin=[0];
fx=[0];
imin=[0];
imaxmin=[0];
ix=[0];
if pc==1
    C=1;
else
    C=3;
end
rf=1;           %proportionality constant of facilitation
fk=-0.01/C;
ri=.14/C^2;     %proportionality constant of inhibition
ik=-0.055/C^2;

%%%%% STAGE I & II %%%%%
for n=1:burst
    %%% STAGE I: Calcium accumulation
    CAmax=C*(1-exp(n*kt))/(1-exp(kt));
    CAmin=C*(1-exp(n*kt))/(1-exp(kt))*exp(kt);
    CAmaxmin=[CAmaxmin CAmin];  % Only CAmin is drawn on the picture
    CAx=[CAx .2*n];
    
    %%% STAGE II: facilitatory/inhibitory substances
    
    %Facilitation
    if n==1
        fmin=CAmin;
    else
        fmin=(rf*(CAmin-CAmaxmin(n-1))+fmin)*exp(fk);
    end
    fmaxmin=[fmaxmin fmin];
    fx=[fx .2*n];
    
    %Inhibition
    imin=(ri*CAmin+imin)*exp(ik);
    imaxmin=[imaxmin imin];
    ix=[ix .2*n];
end
fvalue=fmaxmin;
fxaxis=fx;
ivalue=imaxmin;
ixaxis=ix;
CAvalue=CAmaxmin;
CAxaxis=CAx;

%%% Decay of substances during gaps

base=[.1];
base_gap = 0.1:0.1:gap;
if train > 1
    for de=0.2:0.1:gap-0.1
        base=[base de];
    end
    for cycle=1:train-1
        bfmin=fvalue(length(fvalue))*exp(-bk*base);
        bimin=ivalue(length(ivalue))*exp(-bk*base);
        bcamin=max(CAvalue)*exp(-base); %arbitrary decay number choosen
        fvalue=[fvalue bfmin bfmin(length(bfmin))+fmaxmin];
        ivalue=[ivalue bimin bimin(length(bimin))+imaxmin];
        CAvalue=[CAvalue bcamin bcamin(length(bcamin))+CAmaxmin];
        fxaxis=[fxaxis (gap+.2*burst)*cycle-gap+base (gap+.2*burst)*cycle+fx];
        ixaxis=[ixaxis (gap+.2*burst)*cycle-gap+base (gap+.2*burst)*cycle+ix];
        CAxaxis=[CAxaxis (gap+.2*burst)*cycle-gap+base (gap+.2*burst)*cycle+CAx];
    end
end


%%%%%% Stage III: After Effect %%%%%

%peak effect
pulseno=burst*train;
h=3.6;
peakINH=round(1800*pulseno^h/(130^h+pulseno^h));
peakEXC=round(peakINH/5);
FC=fvalue(length(fvalue)); %final value of facilitation
IC=ivalue(length(ivalue)); %final value of inhibition


% FBC evaluation
if fbc==0
    net=FC-IC;
else
    net=FC;
end

%parameters of onset curve
h1f=3;
h1i=2.5;
t50oEXC=peakEXC*.25;
t50oINH=peakINH*.2;

%parameters of decay curve
h2f=2;
h2i=4;
t50dEXC=peakEXC*2;
t50dINH=peakINH*1.1;

%after effect
aevalue=[];
time=[];
if or (IC >= FC, pc==0)
    for ae=0:3600
        if ae<=peakINH
            I=net*ae^h1i/(t50oINH^h1i+ae^h1i);
            Ip=I;
        else
            I=Ip*t50dINH^h2i/(t50dINH^h2i+(ae-peakINH)^h2i);
        end
        aevalue=[aevalue I];
        time=[time ae];
    end
else
    for ae=0:3600
        if ae<=peakEXC
            F=net*ae^h1f/(t50oEXC^h1f+ae^h1f);
            Fp=F;
        else
            F=Fp*t50dEXC^h2f/(t50dEXC^h2f+(ae-peakEXC)^h2f);
        end
        aevalue=[aevalue F];
        time=[time ae];
    end
end


MAF = HuangModel_Old(time,burst,train,0.16,gap,pc,fbc);


% figure
% subplot(311)
% plot(CAxaxis,CAvalue);
% legend('Ca^{2+}')
% %axis ([0, 1000, 0, 0.5]);
% subplot(312);
% plot(fxaxis,fvalue,'r--');
% hold;
% plot(ixaxis,ivalue);
% legend('Facilitation','Inhibition')
% %axis ([0, 1000, 0, 1.5]);
% hold off;
% subplot(313)
% plot(time,aevalue);
% legend('After Effect')
% %axis ([0, 1200, 0, 15]);


figure
hold on
plot(time,aevalue)
plot(time,MAF)
legend('old_original','old_my')
hold off




