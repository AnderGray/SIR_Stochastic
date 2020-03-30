%%%
%   Quick ecdf, taken from the NasaUQ2020 challenge.
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Liverpool Nasa Lads
%                          Email: ander.gray@liverpool.ac.uk
%
%%%
function [ps,xs] = ecdfNasa(x)
    
    x = rmmissing(x);
    xs = sort(x);
    xs = [xs(1);xs];
    ps = linspace(0,1,length(xs));
    
end