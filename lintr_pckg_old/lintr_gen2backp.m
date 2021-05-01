        function [ys2,var_ep2,as2_mean] = lintr_gen2backp(ys,as,m,n,k,var_ep)
%
%        converts general observations A*X + ep to standard spiked model
%        MX + ep2 with colored noise by backprojecting; no normalization
%        or whitening is transformed
%
%        returns the transformed data, and the normalization matrix and 
%        the variances of the new noise term ep2
%
        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);
%
%        make back-projected data ys2
%
        ys2 = conj(as) .* ys;
        var_ep2 = as2_mean .* var_ep;



        end
%
%
%
%
%
