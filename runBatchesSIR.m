%%%
%
%       Script for running the stochastic SIR model. You may give the simulation the number of
%       batches you wish to run. Will return Ntimes x Nbatches.
%
%%%

addpath("./plottingTools/")
addpath("./Model&Tools/")

Npop = 10000; InInitial = 1;

aplha = 0.2; beta = 1;
V = 100; 

Tsart = 0; Tend = 100;

Nbatches = 200;

Times = linspace(0,5,5000);

[outSn, outIn, outRn] = SIRmc(Tsart, Tend, V, aplha, beta, Npop, InInitial,Nbatches,Times);

outSn = outSn./Npop;
outIn = outIn./Npop;
outRn = outRn./Npop;


figure
set(gcf, 'Position',  [500, 1000, 1000, 800])
p1 = plot(Times,outSn, 'g');
hold on
p2 = plot(Times,outIn, 'r');
p3 = plot(Times,outRn,  'b');

legend('Susceptable','Infected', 'Recovered')
xlim([0 5])

for i=1:length(p1)
    p1(i).Color(4) = 0.6;
    p2(i).Color(4) = 0.6;
    p3(i).Color(4) = 0.6;
end


title("Stochastic SIR model")
xlabel("Time [arb]")
ylabel("% of Population")
set(gca,'FontName','Arial','FontSize',22);


sliceTime(0.7,outSn,outIn,outRn,Times)
