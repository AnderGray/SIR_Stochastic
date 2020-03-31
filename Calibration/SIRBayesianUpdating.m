%%%
%   Script for running Bayesian updating of the SIR stochastic model with Tmcmc
%
%
%                                          Ander Gray
%                                          ander.gray@liverpool.ac.uk
%%% 

%% Location of the model and Extract Data


addpath("../plottingTools/")
addpath("../Model&Tools/")
addpath("./CalibrationData/")

LoadUKData;

Npar = 40;

Nsamples = 120;                 % Samples of the posterior
NBatches = 100;                 % Samples of the stochastic model (Inner loop)

Times = 0:TotalDays;            % Data is given in days. As for prediction at these days

Npop = Population;              % Actual population
Simpop = 10^5;                  % Simulated population

INInitial = TotalDailyInfected(1);


%% Construct prior

alphaBounds = [0.08, 0.3];
betaBounds  = [0.0007, 0.02];
%VBounds     = [50, 200];
V = 1;

%priorV     = @(x) unifpdf(x,VBounds(1),     VBounds(2));
priorAlpha = @(x) unifpdf(x,alphaBounds(1), alphaBounds(2));
priorBeta  = @(x) unifpdf(x,betaBounds(1),  betaBounds(2));


%priorPdf = @(x) priorV(x(1)) * priorAlpha(x(2)) * priorBeta(x(3)) ;  % Uniform independent prior
%priorPdf = @(x) priorAlpha(x(2)) * priorBeta(x(3)) ;  % Uniform independent prior
priorPdf = @(x) priorAlpha(x(1)) * priorBeta(x(2)) ;  % Uniform independent prior
% sampleAlpha = @(N) unifrnd(alphaBounds(1), alphaBounds(2), N,1);
% sampleBeta  = @(N) unifrnd(betaBounds(1),  betaBounds(2), N,1);
% sampleV     = @(N) unifrnd(VBounds(1),     VBounds(2), N,1);
% 
% samplePrior2 = @(N) [sampleV(N), sampleAlpha(N), sampleBeta(N)];

%samplePrior = @(N) LatinPrior(N, [VBounds(1), alphaBounds(1), betaBounds(1)], [VBounds(2), alphaBounds(2), betaBounds(2)]);
samplePrior = @(N) LatinPrior(N, [alphaBounds(1), betaBounds(1)], [alphaBounds(2), betaBounds(2)]);

%% Model Updating
%Model = @(theta) SIRmc(0,TotalDays, theta(1), theta(2), theta(3), Simpop, INInitial, NBatches, Times,Npop);
Model = @(theta) SIRmc(0,TotalDays, V, theta(1), theta(2), Simpop, INInitial, NBatches, Times,Npop);

LogLike = @(theta) LikelihoodArea(Model, theta, TotalDailyInfected, TotalDailyRecovered);

parpool(Npar)

fprintf('Nsamples TMCMC = %d\n', Nsamples);
tic;
SamplesFromPosterior = tmcmc(LogLike, priorPdf, samplePrior, Nsamples, Npar);
toc;

save('SampsPostSIRTmcmc','SamplesFromPosterior');


% Nbatches = Npar;
% NperBatch = Nsamples/Nbatches;
% thetajIn = reshape(SamplesFromPosterior,NperBatch,Nbatches, size(SamplesFromPosterior,2));
% outIn = zeros(5001,Ninner,NperBatch, Nbatches);
% 
% 
% fprintf('Evaluating the model in %d batches\n', Nbatches);
% parfor jk = 1:Nbatches
%     for kk = 1:NperBatch
%        outIn(:,:,kk,jk) = Model(reshape(thetajIn(7kk,jk,:),1,24));
%     end
% end
% 
% simsOut = reshape(outIn,5001,Ninner,Nsamples);
% 
% pboxT = computeBounds(simsOut);
% 
% save('SampsFromPostTmcmc','SamplesFromPosterior', 'pboxT','Ninner','Nsamples');
% 
% pboxF = computeBounds(abs(fft(simsOut)));
% 
% save('SampsFromPostTmcmc','SamplesFromPosterior', 'pboxT','pboxF','Ninner','Nsamples');

function samps = sobolPrior(N, boundsLower,boundsUpper)

    dims = length(boundsLower);
    p = sobolset(dims);
    sobolSamps = net(p,N);

    samps = sobolSamps .* (boundsUpper - boundsLower) + boundsLower;

end

function samps = LatinPrior(N,boundsLower,boundsUpper)

    dims = length(boundsLower);
    LatinSamps = lhsdesign(N,dims);
    samps = LatinSamps .* (boundsUpper - boundsLower) + boundsLower;
end
