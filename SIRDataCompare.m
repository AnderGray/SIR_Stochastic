%%%
%
%       Script for running the stochastic SIR model. You may give the simulation the number of
%       batches you wish to run. Will return Ntimes x Nbatches.
%
%%%

addpath("./plottingTools/")
addpath("./Model&Tools/")
addpath("Calibration/CalibrationData/")

LoadUKData;

Npop = Population;

Simpop = 10^6; InInitial = 1;

aplha = 0.12; beta =  0.005;
V = 1; 

Tsart = 0; Tend = TotalDays+10;

Nbatches = 10;

Times = linspace(Tsart,Tend, 5000);

BinWidth = 0.5;

tic;
%[outSn, outIn, outRn] = SIRmcImportance(Tsart, Tend, V, aplha, beta, Simpop, InInitial, Nbatches,Times, Npop, BinWidth);
[outSn, outIn, outRn] = SIRmc(Tsart, Tend, V, aplha, beta, Simpop, InInitial, Nbatches,Times, Npop);
toc;

figure
set(gcf, 'Position',  [500, 1000, 1000, 800])
p1 = plot(Times,outSn, 'g');
hold on
p2 = plot(Times,outIn, 'r');
p3 = plot(Times,outRn,  'b');

legend('Susceptable','Infected', 'Recovered')

for i=1:length(p1)
    p1(i).Color(4) = 0.6;
    p2(i).Color(4) = 0.6;
    p3(i).Color(4) = 0.6;
end

xlim([0 Tend])
ylim([0 5*10^4])


title("Stochastic SIR model")
xlabel("Time [arb]")
ylabel("Population Numbers")
set(gca,'FontName','Arial','FontSize',22);

plot(1:1:TotalDays,TotalDailyInfected, '*')
plot(1:1:TotalDays, TotalDailyRecovered, '*')

%sliceTime(0.7,outSn,outIn,outRn,Times)
