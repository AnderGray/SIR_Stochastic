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
%       sampled from equivalent distritibutions.  For instance the particle
%       n will undergoes a transition at time t^*. Then at each time step
%       (only used to collect the data) the new transition probability is
%       computed and the transition time t* adjusted.
%
%       TODO: I will write a proper explanation later.
%
%      Edoardo Patelli (edoardo.patelli@strath.ac.uk)%
%%%


%%%
%       Input parameter description
%
%       Vsteps:     Vector of time step using to collect the data. 
%       Mrates:     Matrix of exchange rates 
%                   Mrates(2,1) exchange rate between particle type 1 and 2
%                   Mrates(1,2) exchange rate between particle type 2 and 1....
%                   It is an matrix of halde functions for instance:
%                   
%                   @Vpopulation(alpha * Vpopulation(1) * Vpopulation(2)/ (V * Simpop))
%       Vinitial:   Vector of initial population for particles 
%
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       Output description
%
%       Mout:      output population. Size == length(Vsteps) x NparticleType
%
%
%%%

function [Mout] = MCnonlinear(Vsteps, Mrates, Vinitial)

% NBatches are the number of particles
Nsteps=length(Vsteps); % Number of steps
Ntypes=size(Mrates,1);  % Number of particle type
Nparticles=sum(Vinitial); % total number of particle to simulate

%% Initialisation
% Initialize number of particles
VcurrentPopulation=Vinitial;

Mout=zeros(Nsteps,Nparticles);

% Compute/Update transition probabilities
MratesNow=Mrates(VcurrentPopulation);

TotalRate=cumsum(MratesNow,1); % Probality of transition (equivalent to total cross section)
% Sample transition times

Vtypes=zeros(Nparticles,1);

pStart=1;
for n=1:Ntypes
    pEnd=pStart+Vinitial(n);
    Vtypes(pCount:pEnd)=repmat(n,Vinitial(n),1);
    pStart=pEnd;
    
    VtransitionTime=exprnd(TotalRate(Vtype,end));
end


for n=1:Nsteps  % Loop over the time steps
    TotalRateOld=TotalRate;
    % TODO: Remove the loop and use only vectors and indices
    [Vpos, Vindex]=find(VtransitionTime<Vstep(n));
    
    % Sample reaction type.
    Vsamples=rand(length(Vpos),1);
    

    for p=1:length(Vpos) % loop over the particle type
            % Particle is undergoing a transition in the current time step
            newType=find(TotalRate(Vpos(p),:)>Vsamples(p),'last');
            
            currentType=Vtypes(Vpos(p));
            VcurrentPopulation(currentType)=VcurrentPopulation(currentType)-1;
            
            % Update particle type
            Vtypes(Vpos(p))=newType;
            VcurrentPopulation(newType)=VcurrentPopulation(newType)+1;
  
    end

    
    % Update Transitoon rates
    MratesNow=Mrates(VcurrentPopulation);
    TotalRate=cumsum(MratesNow,1); % Probality of transition (equivalent to total cross section)
    
    % Update transition times
    VtransitionTime(Vpos)=exprnd(TotalRate(Vpos,end));
    
    % Update transition time for non interacting particle
    
    %% This is the only formula that needd to be checked. 
    VtransitionTime(~Vindex)=Vstep(n)+(VtransitionTime(~Vindex)-Vstep(n))* TotalRateOld(~Vindex,end)/TotalRate(~Vindex,end);
    
    % Update counter for the output
    Mout(n,:)=VcurrentPopulation;
    
    % TODO: Instead to add the number of particles in each time step we
    % should add the fraction of time a particle was in the time step. 
    % 
    % e.g. if particle A at time t1 becames B we should add to the counter
    % (Vstep(n-1)+t1)/Vstep(n) and B*(Vstep(n)-t1)/Vstep(n)
    
end

