%%%
%   Produces a plot of the 2nd order distributions of N_R, N_I and N_S from a
%   Monte Carlo simulation of the stochastic SIR model.
%
%   You may also input a single batched simulation and it will plot the ecdf
%
%
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%


function [fig] = sliceTime(t,outSn, outIn, outRn, outT)
   
    Index = find(t >= outT);
    Index = Index(end);
    
    numMc = size(outSn,3);
    
    meanSn = nanmean(outSn(Index,:,:),1);
    meanIn = nanmean(outIn(Index,:,:),1);
    meanRn = nanmean(outRn(Index,:,:),1);
    
    meanSn = meanSn(:);
    meanIn = meanIn(:);
    meanRn = meanRn(:);
    
    f1 = figure;
    set(gcf, 'Position',  [500, 1000, 1000, 800])
    hold on;
    for i = 1:numMc
        [A,B] = ecdf(outSn(Index,:,i));
        p2 = stairs(B,A, '-g');
        title(["t = "; num2str(t)])
        xlabel('Suseptable(t)');ylabel('cdf');set(gca,'FontName','Arial','FontSize',16);
    end
    [A,B] = ecdf(meanSn);
    p2 = stairs(B,A, 'k','LineWidth',4);
    
    f2 = figure;
    set(gcf, 'Position',  [500, 1000, 1000, 800])
    hold on
    for i= 1:numMc        
        [A,B] = ecdf(outIn(Index,:,i));
        p2 = stairs(B,A, '-r');
        title(["t = "; num2str(t)])
        xlabel('Infections(t)');ylabel('cdf');set(gca,'FontName','Arial','FontSize',16);
    end
    [A,B] = ecdf(meanIn);
    p2 = stairs(B,A, 'k','LineWidth',4);
    
    f3 = figure;
    set(gcf, 'Position',  [500, 1000, 1000, 800])
    hold on
    for i = 1:numMc
        [A,B] = ecdf(outRn(Index,:,i));
        p2 = stairs(B,A, '-b');
        title(["t = "; num2str(t)])
        xlabel('Recovered(t)');ylabel('cdf');set(gca,'FontName','Arial','FontSize',16);
    end
    [A,B] = ecdf(meanRn);
    p2 = stairs(B,A, 'k','LineWidth',4);
    
    
    
    
end