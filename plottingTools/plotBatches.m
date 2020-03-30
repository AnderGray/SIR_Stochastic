%%%
%   Plots all batched processes from a single simulation
%
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%   t:      Start time of time lapse
%   Pred:   2nd order distribution. Dims == times x Inners x Outters
%   times:  Vector of time points
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%%%
function [fig] = plotBatches(outSn, outIn, outRn, outT, PlotAlpha)
    
    figure
    set(gcf, 'Position',  [500, 1000, 1000, 800])
  
   
    
    p1 = plot(outT,outSn, 'g');
    legend('Susceptable')
    hold on
    p2 = plot(outT,outIn, 'r');
    legend('Infected')
    p3 = plot(outT,outRn,  'b');
    legend('Recovered')

    %legend('Susceptable','Infected', 'Recovered')
    xlim([0 15])

    for i=1:length(p1)
        p1(i).Color(4) = PlotAlpha;
        p2(i).Color(4) = PlotAlpha;
        p3(i).Color(4) = PlotAlpha;
    end


    title("Stochastic SIR model")
    xlabel("Time [arb]")
    ylabel("% of Population")
    set(gca,'FontName','Arial','FontSize',22);
    
    fig = gca;
end