%%%
%   A simple Susceptable - Infected - Recovered stochatsic model for
%   epidemic modelling using non linear Monte Carlo Simulation
%
%       This method derives from the Monte Carlo method for particle transport.
%       See for more details:
%       https://gnssn.iaea.org/NSNI/Shared%20Documents/OPEN%20Shared%20Files/MonteCarloParticleTransportMethodsNeutronAndPhotonCalculations.pdf
%
%       It allows to simulate the fate of particles (agent) with non-linear
%       transition rates (i.e. where the probabilities of interaction
%       depends on the number of particle present).
%
%       See references:
%       * https://doi.org/10.1016/S0306-4549(00)00082-7
%       * https://doi.org/10.1016/S0306-4549(02)00072-5
%       * https://doi.org/10.1016/j.anucene.2006.11.011
%
%       There are two ways to perform Monte Carlo particle transport with
%       non-linear parameters:
%       1) Perform the analysis on small time steps where the rates are
%           assumed to be constant and hence sample if the particle is
%           undergoing a transition in the time step or not. Then update
%           the parameters using the number of particles present and
%           process the next time step.
%           This approach is simple and it does not require to change the
%           linear Monte Carlo particle transport approach. Only adding a
%           loop over the time steps.
%
%       2) In the second approach the transition time of each particle is
%       samples from equivalent distritibutions.  For instance the particle
%       n will undergoes a transition at time t^*. Then at each time step
%       (only used to collect the data) the new transition probability is
%       computed and the transition time t* adjusted.
%
%       TODO: I will write a proper explanation later.
%
%      Code adapted from the the SIRmc model from Ander Gray.
%
%      Edoardo Patelli (edoardo.patelli@strath.ac.uk)%
%%%


%%%
%       Input parameter description
%
%       Tstart:     Initial time of simultaion
%       Tend:       Final time of simulation
%       V:          Spacial parameter. Models social distancing
%       alpha:      Infection rate parameter
%       beta:       Recovery rate parameter
%       Simpop:     Total simulated population
%       INInitial:  Initial number of infected population
%       Nbatches:   Number of times the stochastic model will be run
%       outTpoints: Points in time to find predictions. Will interpolate at these points.
%       Npop:       True population
%
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       Output description
%
%       outSn:      output Susceptable population. Size == length(outTpoints) x Nbatches
%       outIn:      output Infected population.    Size == length(outTpoints) x Nbatches
%       outRn:      output Recovered population.   Size == length(outTpoints) x Nbatches
%
%
%%%

function [outSn, outIn, outRn] = SIRmcNL(Tstart, Tend, V, alpha, beta, Simpop, INInitial, Nbatches, outTpoints, Npop)

% NBatches are the number of particles
Pfraction = 1/100;

NumTpoints = length(outTpoints);

outSn = zeros(NumTpoints, Nbatches);
outIn = zeros(NumTpoints, Nbatches);
outRn = zeros(NumTpoints, Nbatches);

for i =1:Nbatches
    N_I = INInitial;
    N_S = Simpop - N_I;
    N_R = 0;
    
    outSnBatch = [N_S]; outInBatch = [N_I];
    outRnBatch = [N_R]; outTsBatch = [Tstart];
    
    t = Tstart;
    
    while t < Tend
        if N_I == 0; break; end
        
        w1 = alpha * N_S * N_I / (V * Simpop);         %   Rate at which infections occur
        w2 = beta * N_I;                    %   Rate at which cures occur
        
        w = w1 + w2;                        %   Rate at which reactions occur
        
        % Time step selected to have on average Pfraction (e.g. 1/100) events in dt.
        
        Tstep=expinv(Pfraction,w);          %   Calculate time step
        Tstar=exprnd(w);                    %   Samples of the inverse cdf of an exp(w)
        
        t = t + Tstep;                      %   Time step
        
        if Tstar>Tstep
            % In this case nothing happen. The transition is after the
            % considered time step.
        else
            % The particle undergoes a transition
            
            if rand() < w1/w                    %   Which reaction?
                
                N_S = N_S - 1;                  %   Infection
                N_I = N_I + 1;
                
            else
                
                N_R = N_R + 1;                  %   Recovery
                N_I = N_I - 1;
                
            end
        end
        
        % Not a very efficient way to store the results
        outTsBatch = [outTsBatch; t];
        outSnBatch = [outSnBatch; N_S];
        outInBatch = [outInBatch; N_I];
        outRnBatch = [outRnBatch; N_R];
        
    end
    
    outSn(:,i) = interp1(outTsBatch,outSnBatch,outTpoints, 'nearest');
    outRn(:,i) = interp1(outTsBatch,outRnBatch,outTpoints, 'nearest');
    outIn(:,i) = interp1(outTsBatch,outInBatch,outTpoints, 'nearest');
    
end

outSn = outSn/Simpop * Npop;
outIn = outIn/Simpop * Npop;
outRn = outRn/Simpop * Npop;

end