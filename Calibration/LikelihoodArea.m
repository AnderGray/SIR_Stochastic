function LogL = LikelihoodArea(model, theta, DInfected, DRecovered)

    LogLInf = zeros(size(DInfected,1),1);
    LogLRec = zeros(size(DRecovered,1),1);
    [~, predInf, predRec] = model(theta);
    epsilon = 0.1;  % "width factor" of the likeihood
                    % Should be [10^-3 , 10^-1] according to Metric paper.
                    % The smaller it is the more spiky the posterior
                    % but the more compute we need to sample
    
    predInf(isnan(predInf))  =   0;
    predRec(isnan(predRec))  =   max(max(predRec));
                    
    for i = 1:size(DInfected,1)
        
        stochasticDistanceInf = areaMe(predInf(i,:),DInfected(i,:));
        stochasticDistanceRec = areaMe(predRec(i,:),DRecovered(i,:));
        
        LogLInf(i) = -(1/epsilon)^2 * stochasticDistanceInf^2;
        LogLRec(i) = -(1/epsilon)^2 * stochasticDistanceRec^2;
        
    end
    LogL = sum(LogLInf) + sum(LogLRec);
end
