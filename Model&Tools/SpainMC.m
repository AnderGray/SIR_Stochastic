%%%
%   A simple Susceptable - Infected - Recovered stochatsic model for epidemic modelling
%
%       This method has sometimes been called Dynmamic Monte Carlo, and is 
%       used in reaction chemistry for predicting the population of some chemical
%       compounds at a future time, for some copuled reaction pathways with some rate
%       at which the reactions occur; and that the chemicals are well mixed.
%
%       The problem may also be solve with coupled ordinary differential equations, 
%       but becomes difficult if the ODE's are too stiff. Dynamic MC does not have
%       this problem.
%
%       The problem is almost idential for an SIR model, and works in the following way.
%       The initial populations of N_I, N_S are given, with N_R = 0; and also
%       the rates at which people are infected and recovered (InfRate, CureRate).
%       Starting from t=0, either an infection or a recovery will occor. The time to 
%       one of these reactions is sampled from an exponential distribution, whose parameter
%       is calculated from N_I, N_S, InfRate and CureRate. The time is stepped by the sample
%       and then one of the reactions is sampled. The populations is then ajusted according
%       to which reaction has occured. This is repeated until t = Tfinal or until N_I == 0.
%
%       Simple!
%           
%       Because this is a stochastic simulation, the entire simulation may be re-run (batched)
%       a number of times calculating the variation due to the stochasticity.
%       In this model this is the Nbatches parameter.
%
%       We can also model social distancing. This is parameter V, and it reduces
%       the rate of infection.
%
%
%       If you find that this model is too slow, we can run it in parallel, 
%       implement importance sampling, or both.
%
%       For more info: https://en.wikipedia.org/wiki/Gillespie_algorithm#Another_example:_The_SIR_epidemic_without_vital_dynamics
%
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Enrique Miralles & Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
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

function [outSn, outIn, outRn, outQn, outLn, outHn, outUn, outHUn, outFn] = SpainMC(Tstart, Tend, QuarentineTime, rates, InitialConditions, Nbatches, outTpoints, Simpop, Npop)
    
    InitialConditions = round(InitialConditions,0);
    
    Simpop = sum(InitialConditions);
    
    betaBefore   = rates(1);
    betaAfter = rates(2);
    gamma1 = rates(3);
    gamma2 = rates(4);
    gamma3 = rates(5);
    alpha1 = rates(6);
    alpha2 = rates(7);
    alpha3 = rates(8);
    eta    = rates(9);
    d1     = rates(10);
    d2     = rates(11);
    delta  = rates(12);
    tau    = rates(13);
    
    % Tau is a parameter which controls how many people are in quaretine as apposed to being
    % susceptable.
    
    InputedPop = sum(InitialConditions);
    
    if InputedPop ~= Simpop
       disp("There's an error") 
    end
    
    Sinitial = InitialConditions(1);
    Iinitial = InitialConditions(2);
    Rinitial = InitialConditions(3);
    Qinitial = InitialConditions(4);
    Linitial = InitialConditions(5);
    Hinitial = InitialConditions(6);
    Uinitial = InitialConditions(7);
    HUinitial = InitialConditions(8);
    Finitial = InitialConditions(9);

    
    NumTpoints = length(outTpoints);
    
    outSn  = zeros(NumTpoints, Nbatches);
    outIn  = zeros(NumTpoints, Nbatches);
    outRn  = zeros(NumTpoints, Nbatches);
    outQn  = zeros(NumTpoints, Nbatches);
    outLn  = zeros(NumTpoints, Nbatches);
    outHn  = zeros(NumTpoints, Nbatches);
    outUn  = zeros(NumTpoints, Nbatches);
    outHUn = zeros(NumTpoints, Nbatches);
    outFn  = zeros(NumTpoints, Nbatches);
    
    
    for i =1:Nbatches
        
        beta = betaBefore;
        tauIn = 1;
        deltaIn = 0;
        
        outSnBatch = [Sinitial]; outInBatch = [Iinitial]; 
        outRnBatch = [Rinitial]; outQnBatch = [Qinitial]; 
        outLnBatch = [Linitial]; outHnBatch = [Hinitial]; 
        outUnBatch = [Uinitial]; outHUnBatch = [HUinitial]; 
        outFnBatch = [Finitial]; 
        
        N_S = Sinitial; N_I = Iinitial; N_R = Rinitial;
        N_Q = Qinitial; N_L = Linitial; N_H = Hinitial; 
        N_U = Uinitial; N_HU=HUinitial; N_F = Finitial;
        
        outTsBatch = [Tstart];

        t = Tstart;

        while t < Tend
            if N_I == 0; break; end
        
            if t >= QuarentineTime
                beta = betaAfter;
                tauIn= tau;
                deltaIn = delta;
            end

            wSL = beta * N_I/Simpop * N_S;
            wSQ = deltaIn * N_S;                    
            wQS = tauIn * N_Q;
            wLI = gamma1 * N_L;
            wIH = gamma2 * N_I;
            wIR = alpha1 * N_I;
            wHU = gamma3 * N_H;
            wHF = d1     * N_H;
            wHR = alpha2 * N_H;
            wUHU= alpha3 * N_U;
            wUF = d2     * N_U;
            wHUR= eta    * N_HU;

            w = wSL + wSQ + wQS + wLI + wIH + wIR + wHU + wHF + wHR + wUHU + wUF + wHUR; 
            
            dt = -log(rand())/w;              

            t = t + dt;                     
            
            probs = [wSL/w, wSQ/w, wQS/w, wLI/w, wIH/w, wIR/w, wHU/w, wHF/w,...
                wHR/w, wUHU/w, wUF/w, wHUR/w];
            
            cdf = cumsum(probs);
            Index = find(rand() <= cdf);
            Event = Index(1);
            
            if Event == 1

                N_S = N_S - 1;                  %   Infection
                N_L = N_L + 1;

            elseif Event == 2

                N_S = N_S - 1;                 
                N_Q = N_Q + 1;
                
            elseif Event == 3

                N_Q = N_Q - 1;
                N_S = N_S + 1;                 
                
            elseif Event == 4

                N_L = N_L - 1;                 
                N_I = N_I + 1;
                
            elseif Event == 5

                N_I = N_I - 1;                 
                N_H = N_H + 1;
                
            elseif Event == 6

                N_I = N_I - 1;                 
                N_R = N_R + 1;
                
            elseif Event == 7

                N_H = N_H - 1;                 
                N_U = N_U + 1;

            elseif Event == 8

                N_H = N_H - 1;                 
                N_F = N_F + 1;               
                
            elseif Event == 9

                N_H = N_H - 1;                 
                N_R = N_R + 1;
                
            elseif Event == 10

                N_U = N_U - 1;                 
                N_HU = N_HU + 1;

            elseif Event == 11

                N_U = N_U - 1;                 
                N_F = N_F + 1;
                
            elseif Event == 12

                N_HU = N_HU - 1;                 
                N_R = N_R + 1;

            end

            outTsBatch = [outTsBatch; t];
            
            outSnBatch = [outSnBatch; N_S];
            outInBatch = [outInBatch; N_I];
            outRnBatch = [outRnBatch; N_R];
            outQnBatch = [outQnBatch; N_Q];
            outLnBatch = [outLnBatch; N_L];
            outHnBatch = [outHnBatch; N_H];
            outUnBatch = [outUnBatch; N_U];
            outHUnBatch= [outHUnBatch; N_HU];
            outFnBatch = [outFnBatch; N_F];

            
        end
        
        outSn(:,i) = interp1(outTsBatch,outSnBatch,outTpoints, 'nearest');
        outRn(:,i) = interp1(outTsBatch,outRnBatch,outTpoints, 'nearest');
        outIn(:,i) = interp1(outTsBatch,outInBatch,outTpoints, 'nearest');
        outQn(:,i) = interp1(outTsBatch,outQnBatch,outTpoints, 'nearest');
        outLn(:,i) = interp1(outTsBatch,outLnBatch,outTpoints, 'nearest');
        outHn(:,i) = interp1(outTsBatch,outHnBatch,outTpoints, 'nearest');
        outUn(:,i) = interp1(outTsBatch,outUnBatch,outTpoints, 'nearest');
        outHUn(:,i) = interp1(outTsBatch,outHUnBatch,outTpoints, 'nearest');
        outFn(:,i) = interp1(outTsBatch,outFnBatch,outTpoints, 'nearest');
    
    end
    
    outSn = outSn/Simpop * Npop;
    outIn = outIn/Simpop * Npop;
    outRn = outRn/Simpop * Npop;
    outQn = outQn/Simpop * Npop;
    outLn = outLn/Simpop * Npop;
    outHn = outHn/Simpop * Npop;
    outUn = outUn/Simpop * Npop;
    outHUn = outHUn/Simpop * Npop;
    outFn = outFn/Simpop * Npop;
    
    
    
end