%%%
%   This script is an altered version of the vine.m script for
%   generating a random correlation matrix
%
%   The scripts generates a correlation matrix from partial
%   correlations, as described in:
%
%       Joe, Harry, "Generating random correlation matrices based on partial correlations." 
%           Journal of Multivariate Analysis 97.10 (2006): 2177-2189.
%
%   The partial correlations may vary between [-1, 1] and we need
%   nchoosek(d,2) of them. When these parameters are randomised, a
%   random correlation matrix with (nearly) uniformly distributed
%   correlations is produced
%
%   It is the partialCorrs which will be calibrated in the bayesian
%   updating
%
%   D           == dimension of the correlation matrix
%   partialCorr == the partial correlations, nchoosek(d,2) of them
%
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%
function S = generateCorr(d, partialCorrs)
    
    P = zeros(d);           %// storing partial correlations
    S = eye(d);
    ll=1;
    for k = 1:d-1
        for i = k+1:d
            P(k,i) = partialCorrs(ll);
            p = P(k,i);
            for l = (k-1):-1:1 %// converting partial correlation to raw correlation
                p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
            end
            S(k,i) = p;
            S(i,k) = p;
            ll=ll+1;
        end
    end

end