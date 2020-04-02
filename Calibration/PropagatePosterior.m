%%%S
%       Script for running uncertainty propagation through stochastic SIR model.
%
%       The uncertain parameters are the infection rate, recovery rate and spacial
%       parameter.
%
%       Input distributions are Gaussian, but may be anything. Dependency is
%       is modelled using a gaussian copula, which is sampled with Cholesky decomposition.
%       You may define the correlations in terms of the partial correlations, which
%       may take any value in [-1, 1], independently.
%
%       Because we have both aleatory (stochastic model) and epistemic uncertainty (uknown rates),
%       the output of this process will be a 2nd order distribution, a distribution of
%       distributions. We will have one for every point in time. You may take the bounds of the 2nd
%       order distribution and create a p-box. 
%
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%

addpath("./plottingTools/")
addpath("./Model&Tools/")

load("SampsPostSIRTmcmc.mat")


addpath("../plottingTools/")
addpath("../Model&Tools/")
addpath("./CalibrationData/")

LoadUKData;


Npop = Population;              % Actual population
Simpop = 10^5;                  % Simulated population


Npar = 36;              % Number of parpools


Nmc = size(SamplesFromPosterior,1);                 % Samples of the uncertainty(Outter loop)
NBatches = 100;           % Samples of the stochastic model (Inner loop)


%%%% 
% Evaluating the stochastic model for random inputs

numPoints = 5000;

alpha = SamplesFromPosterior(:,1);
beta =  SamplesFromPosterior(:,2);

Tsart = 0; Tend = TotalDays+10;
Times = linspace(Tsart,Tend,numPoints);

V = 1;

InInitial = TotalDailyInfected(1);

outSn = zeros(numPoints,NBatches, Nmc);
outIn = zeros(numPoints,NBatches, Nmc);
outRn = zeros(numPoints,NBatches, Nmc);

parpool(Npar)

tic;
parfor i=1:Nmc
   
    [outSn(:,:,i), outIn(:,:,i), outRn(:,:,i)] = SIRmc(Tsart, Tend, V, alpha(i), beta(i), Simpop, InInitial,NBatches,Times, Npop);
    
end
toc

outT = Times;

meanSn = mean(nanmean(outSn,3),2);
meanIn = mean(nanmean(outIn,3),2);
meanRn = mean(nanmean(outRn,3),2);

ff = figure;
set(gcf, 'Position',  [500, 1000, 1000, 800]);
hold on

% Plot all processes
for i = 1:Nmc
    p1 = plot(outT,outSn(:,:,i), 'g');
    p2 = plot(outT,outIn(:,:,i), 'r');
    p3 = plot(outT,outRn(:,:,i),  'b');
    
    for ii=1:length(p1)
        p1(ii).Color(4) = 0.2;
        p2(ii).Color(4) = 0.2;
        p3(ii).Color(4) = 0.2;
    end
    
end

p1 = plot(outT,meanSn, 'k', 'LineWidth',3);
p2 = plot(outT,meanIn, 'k', 'LineWidth',3);
p3 = plot(outT,meanRn, 'k', 'LineWidth',3);

ylim([0 5*10^4])
xlim([Tsart Tend])

plot(1:1:TotalDays,TotalDailyInfected, '*')
plot(1:1:TotalDays, TotalDailyRecovered, '*')

title("Stochastic SIR model")
xlabel("Time [arb]")
ylabel("Number of people")
set(gca,'FontName','Arial','FontSize',22);

saveas(ff,'Process.png')

% Plot p-box at a slice
sliceTimeSave(61,outSn,outIn,outRn,Times);

save('simOut','outSn','outIn','outRn','Times')

% P-box time lapse
%rollingPbox(0, outIn,Times);

%outInBounds = computeBounds(outIn);
%rollingPboxBounds(0,outInBounds,outIn,Times);


