%%%
%
%       Script for running the stochastic SIR model. You may give the simulation the number of
%       batches you wish to run. Will return Ntimes x Nbatches.
%
%%%

addpath("./plottingTools/")
addpath("./Model&Tools/")

Npop = 47*10^6; 

Simpop = 10^5; 

p1 = 5/100; p2 = 6.5/100; p3 = 4.5/100; p4 = 55/100;

%[betaBefore, betaAfter, gamma1, gamma2, gamma3, alpha1, alpha2, alpha3, eta, d1, d2, delta, tau]
rates = [0.59, 0.59, 1/5.2, p2/5.8, p3/1, (1-p1)/14, (1-p2-p3)/7, (1-p4)/14, 1/6, p2/7.5, p4/8, 70.36/100, 0.1];

Tsart = 0; Tend = 270;
QuarentineT = 20;

TauChangeTimes = [60, 90, 170,250];
TauChangeRates = [0.5, 0.1, 1,1];

Nbatches = 1;

Times = linspace(Tsart,Tend,5000);

In = 15*10^3;
Ln = 11*10^3;
Sn = Npop - Ln - In;
Rn = 0;
Qn = 0;
Hn = 0;
Un = 0;
HUn= 0;
Fn = 0;

InInitial = [Sn, In, Rn, Qn, Ln, Hn, Un, HUn, Fn]/Npop * Simpop; 

[outSn, outIn, outRn, outQn, outLn, outHn, outUn, outHUn, outFn] = ValenciaModel(Tsart, Tend, QuarentineT,TauChangeTimes, rates,TauChangeRates, InInitial, Nbatches,Times, Simpop, Npop);



figure
set(gcf, 'Position',  [500, 1000, 1000, 800])
p1 = plot(Times,outSn, 'g');
hold on
p2 = plot(Times,outIn, 'r');
p3 = plot(Times,outRn, 'b');
p2 = plot(Times,outQn, 'Color',[1,0.5,0.5]);
p3 = plot(Times,outLn, 'y');
p2 = plot(Times,outHn, 'c');
p3 = plot(Times,outUn, 'm');
p2 = plot(Times,outHUn,'Color', [0.5,0.5,1]);
p3 = plot(Times,outFn, 'k');

legend('Susceptable','Infected', 'Recovered', 'Q', 'L', 'H', 'U', 'HU', 'F')
xlim([0 Tend])
% 
% for i=1:length(p1)
%     p1(i).Color(4) = 0.6;
%     p2(i).Color(4) = 0.6;
%     p3(i).Color(4) = 0.6;
% end


title("Stochastic SIR model")
xlabel("Time [arb]")
ylabel("Number of Population")
set(gca,'FontName','Arial','FontSize',22);


%sliceTime(0.7,outSn,outIn,outRn,Times)
