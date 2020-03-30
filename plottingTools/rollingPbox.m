%%%
%   Produces a time lapse of the 2nd order MC propagation, with
%   the mean distribution also plotted.
%
%   Inputs:
%
%   t:      Start time of time lapse
%   Pred:   2nd order distribution. Dims == times x Inners x Outters
%   times:  Vector of time points
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%
function out = rollingPbox(t,Pred,times)
    
    Index = find(t >= times);
    Index = Index(end);
    
    f1 = figure;
    set(gcf, 'Position',  [500, 1000, 1000, 800])
    
    for j =Index:10:length(times)
        
       hold on
       xlabel('y(t)');ylabel('cdf');set(gca,'FontName','Arial','FontSize',16);
        
       for i = 1:size(Pred,3)
            [A,B] = ecdf(Pred(j,:,i)');
            p2 = stairs(B,A, '-b');
            %p1.Color(4) = 0.2;
       end
        
        meanPred = nanmean(Pred(j,:,:),1);
        meanPred = meanPred(:);
        [A,B] = ecdf(meanPred);
        p1 = stairs(B,A, '-r','LineWidth', 4);
        %p1.Color(4) = 0.2;
        title(['t = ',num2str(times(j))])
        
        pause(0.000001)
        c1 = get(f1,'children');
        delete(c1)
    end
end