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


Nmc = 50;                 % Samples of the uncertainty(Outter loop)
NBatches = 100;           % Samples of the stochastic model (Inner loop)

%%% 
%   Generating correlated uncertain inputs

partCor1 = 0; partCor2 = 0; partCor3 = 0;                      % Define partial correlations

CorMatrix = generateCorr(3, [partCor1, partCor2, partCor3]);   % Generates a correlation matrix

a = chol(CorMatrix,'lower');           % Cholesky decomp method for generating
z = randn(Nmc, 3);                     % correlated random normals
x = z * a';

u = normcdf(x);                        % Samples from the Gaussian copula

aplha = norminv(u(:,1), 0.2,0.005);
beta  = norminv(u(:,2), 1,0.005);
V     = norminv(u(:,3),100,0.1);

%
%%%%

%%%% 
% Evaluating the stochastic model for random inputs

numPoints = 5000;

Tsart = 0; Tend = 7;
Times = linspace(Tsart,Tend,numPoints);

Npop = 10000; InInitial = 1;

outSn = zeros(numPoints,NBatches, Nmc);
outIn = zeros(numPoints,NBatches, Nmc);
outRn = zeros(numPoints,NBatches, Nmc);


if ~isempty(gcp('nocreate'))        % Create parpool if not already
    parpool(4)
end

parfor i=1:Nmc
   
    [outSn(:,:,i), outIn(:,:,i), outRn(:,:,i)] = SIRmc(Tsart, Tend, V(i), aplha(i), beta(i), Npop, InInitial,NBatches,Times);
    
end

outT = Times;

meanSn = mean(nanmean(outSn,3),2);
meanIn = mean(nanmean(outIn,3),2);
meanRn = mean(nanmean(outRn,3),2);

figure
set(gcf, 'Position',  [500, 1000, 1000, 800])
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

xlim([Tsart Tend])

title("Stochastic SIR model")
xlabel("Time [arb]")
ylabel("Number of people")
set(gca,'FontName','Arial','FontSize',22);


% Plot p-box at a slice
sliceTime(0.7,outSn,outIn,outRn,Times)

% P-box time lapse
rollingPbox(0, outIn,Times);

%outInBounds = computeBounds(outIn);
%rollingPboxBounds(0,outInBounds,outIn,Times);


