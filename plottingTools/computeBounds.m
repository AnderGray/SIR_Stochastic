%%%
%   Computes p-box bounds from a 2nd order Distribution. Can introduce a 
%   vector of 2nd order distributiouns. Eg, the first dimension of Pred may be
%   a indicie in a stochastic process.
%
%                  Institute for risk and uncertainty, University of Liverpool
%
%                          Author: Ander Gray
%                          Email: ander.gray@liverpool.ac.uk
%
%%%


function Pbox = computeBounds(Pred)
    
    ts = size(Pred,1);
    Ninners = size(Pred,2);
    Nouters = size(Pred,3);
    Pbox = zeros(ts,Ninners,2);

    for i = 1:ts
        [~,LB] = ecdfNasa(Pred(i,:,1)'); 
        UB = LB;
        for j = 2:Nouters
           [~,B1] = ecdfNasa(Pred(i,:,j)'); 
           LB = min(LB,B1);
           UB = max(UB,B1);
        end
        Pbox(i,:,1) = LB(2:end);
        Pbox(i,:,2) = UB(2:end);
    end
end