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
%                          Author: Ander Gray
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

function [outSn, outIn, outRn] = SIRmcImportance(Tstart, Tend, V, alpha, beta, Simpop, INInitial, Nbatches, outTpoints, Npop, BinWidth)
    
    
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
        
        weightBinWidth = BinWidth;
        
        while t < Tend
            if N_I <= 0; break; end

            w1 = alpha * N_S * N_I / (V * Simpop);         %   Rate at which infections occur
            w2 = beta * N_I;                    %   Rate at which cures occur
            
            w = w1 + w2;                        %   Rate at which reactions occur
            
            probInf = w1/w;
            probRec = w2/w;
            
            %probInf = w1;
            %probRec = w2;
            
            % Missing some steps here... 
            
            % iii) Construct weights from a weighting functions (do we define?). Must be integers
            % iv)  Construct weighted probabilities  pwi = pi/wi
            % v)   Ptw = sum(pwi)
            % vi)  sample events with pwi. Make reaction i happen wi times
            % vii) sample timestep with Ptw
            
            wInf = ceil(probInf/weightBinWidth);
            wRec = ceil(probRec/weightBinWidth);
            
            wInf = min(wInf, outSnBatch(end));
            %wRec = 1;
            
            
            %probInfImport = (probInf/wInf);
            %probRecImport = (probRec/wRec);
            
            if probInf == 0
                wInfImport = 0;
                wRecImport = 1;
            elseif probRec == 0
                wInfImport = 1;
                wRecImport = 0;
            else
                wInfImport = (probInf/wInf) * w;
                wRecImport = (probRec/wRec) * w;
            end
            % This is wrong
            %PTweight = probInfImport + probRecImport;
            
            wPTweight = wInfImport + wRecImport;
            
            if  rand() <  wInfImport/wPTweight                    %   Which reaction?
            
                for jj = 1:wInf
                    N_S = N_S - 1;                  %   Infection
                    N_I = N_I + 1;
                end

            else
                for jj = 1:wRec
                    N_R = N_R + 1;                  %   Recovery
                    N_I = N_I - 1;
                end
            end
            
            
            %dt = -log(rand())/PTweight;         %   Samples of the inverse cdf of an exp(w)
            dt = -log(rand())/wPTweight;         %   Samples of the inverse cdf of an exp(w)

            t = t + dt;                         %   Time step

            outTsBatch = [outTsBatch; t];
            outSnBatch = [outSnBatch; N_S];
            outInBatch = [outInBatch; N_I];
            outRnBatch = [outRnBatch; N_R];

        end
        
        outSn(:,i) = interp1(outTsBatch,outSnBatch,outTpoints, 'nearest','extrap');
        outRn(:,i) = interp1(outTsBatch,outRnBatch,outTpoints, 'nearest','extrap');
        outIn(:,i) = interp1(outTsBatch,outInBatch,outTpoints, 'nearest','extrap');
    
    end
    
    outSn = outSn/Simpop * Npop;
    outIn = outIn/Simpop * Npop;
    outRn = outRn/Simpop * Npop;
    
end