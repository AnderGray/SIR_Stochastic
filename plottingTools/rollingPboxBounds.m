%%%
%   Produces a time lapse of the bounds of the 2nd order MC propagation, with
%   the mean distribution also plotted.
%   
%   Note: This doesn't quite work due to the NaN's
%
%   Inputs:
%
%   t:            Start time of time lapse
%   PredBounds:   p-box bounds. Dims == times x Inners x 2
%   Pred:         2nd order distribution. Dims == times x Inners x Outters
%   times:        Vector of time points
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%
function out = rollingPboxBounds(t,PredBounds, Pred,times)
    
    Index = find(t >= times);
    Index = Index(end);
    
    f1 = figure;
    set(gcf, 'Position',  [500, 1000, 1000, 800])
    xlabel('y(t)');ylabel('cdf');set(gca,'FontName','Arial','FontSize',16);
    %prompt = 'What is the original value? ';
    %x = input(prompt);
    for j =Index:3:length(times)
       hold on
       
        
 
        [A1,B1] = ecdfNasa(PredBounds(j,:,1)');
        [A2,B2] = ecdfNasa(PredBounds(j,:,2)');
        p2 = stairs(B1,A1, '-b');
        p2 = stairs(B2,A2, '-b');
        %p1.Color(4) = 0.2;
        
        
        A1 = A1(2:end);
        A2 = A2(1:end-1);

        B1 = B1(2:end);
        B2 = B2(2:end);

        Xs = sort([B1;B1;B2;B2]);
        nums = length(Xs);
        Ylb = zeros(nums,1);
        Yub = zeros(nums,1);

        for i = 1:2:nums

            indLB = find(Xs(i) <= B1);
            indUB = find(Xs(i) >= B2);

            if ~isempty(indLB)
                Ylb(i) = A2(indLB(1)); 
                if ismember(Xs(i),B1) 
                    Ylb(i+1) = A1(indLB(1));
                else
                    Ylb(i+1) = A2(indLB(1));
                end
            else 
                Ylb(i) = 1; 
                Ylb(i+1)=1;
            end

            if ~isempty(indUB)
                Yub(i+1) = A1(indUB(end)); 
                if ismember(Xs(i), B1) 
                    Yub(i) = A1(indUB(end));
                else
                    Yub(i) = A2(indUB(end));
                end
            else
                Yub(i) = 0; 
                Yub(i+1)=0;
            end

        end
        cols = [0.9, 0.7, 0.3];

        opts={'EdgeColor', 'none',...
              'FaceColor', cols};
        fh = fill_between(Xs,Ylb,Yub,[],opts{:});
        
    
        %xlim([-0.1 0.1]);
        %ylim([0 1]);
        
            xlim([min(PredBounds(:)) max(PredBounds(:))]);
            ylim([0 1]);


        
%         Xs = sort([B1;B2]);
%         nums = length(Xs);
%         Ylb = zeros(nums,1);
%         Yub = zeros(nums,1);
% 
%         for i = 1:nums
% 
%             indLB = find(Xs(i) <= B1);
%             indUB = find(Xs(i) >= B2);
% 
%             if ~isempty(indLB); Ylb(i) = A2(indLB(1));else; Ylb(i) = 1;end
%             if ~isempty(indUB); Yub(i) = A1(indUB(end));else; Yub(i) = 0;end
% 
%         end
%         
%         
%         fill_color = [.929 .694 .125*5];
%         
%         opts={'EdgeColor', 'none',...
%         'FaceColor', fill_color};
%         fh = fill_between(Xs,Ylb,Yub,[],opts{:});         

        
        meanPred = nanmean(Pred(j,:,:),1);
        meanPred = meanPred(:);
        % Plot exp
        [A,B] = ecdf(meanPred);
        p1 = stairs(B,A, '-r', 'LineWidth',3);
        %p1.Color(4) = 0.2;
        title(['t = ',num2str(times(j))])
        
        pause(0.001)
        c1 = get(f1,'children');
        delete(c1)
    end
end