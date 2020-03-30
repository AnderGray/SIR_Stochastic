%%%
%
%       Script for running the stochastic SIR model
%       Plots the batches from a single stochastic simulation, as they are
%       calculated
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%%%


figure
set(gcf, 'Position',  [500, 1000, 1000, 800])

for i= 1:100
    
    Npop = 1000; InInitial = 1;

    aplha = 2; beta = 1;

    V = 100; 
    outT = linspace(0,10,1000);
    Tsart = 0; Tend = 100;


    [outSn, outIn, outRn] = SIRmc(Tsart, Tend, V, aplha, beta, Npop, InInitial,1,outT);
    
    outSn = outSn./Npop;
    outIn = outIn./Npop;
    outRn = outRn./Npop;
    

    p1 = plot(outT,outSn, 'g');
    p1.Color(4) = 0.2;
    hold on
    p2 = plot(outT,outIn, 'r');
    p2.Color(4) = 0.2;
    p3 = plot(outT,outRn,  'b');
    p3.Color(4) = 0.2;
    xlim([0 15])
    
    
    title("Stochastic SIR model")
    xlabel("Time [arb]")
    ylabel("% of Population")
    legend('Susceptable','Infected', 'Recovered')
    set(gca,'FontName','Arial','FontSize',22);
    drawnow();
    pause(0.01)
end