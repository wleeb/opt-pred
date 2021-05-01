        function [ys4,wvals] = lintr_gen2whi2(ys,as,m,n,k,var_ep,as2_mean)
%
%        converts general observations A*X + ep to standard spiked model
%        with whitened noise by backprojecting, normalizing, and whitening
%
%        returns the whitened data, and the diagonal whitening transformation
%
%
%        make back-projected data ys2, normalized backprojected ys3,
%        and effective noise variances for each one
%
        ys2 = conj(as) .* ys;
        ys3 = ys2 ./ repmat(as2_mean,1,n);
%
        var_ep2 = as2_mean .* var_ep;
        var_ep3 = var_ep2 ./ as2_mean.^2;

%
%        now whiten the effective noise
%
        wvals = 1./sqrt(var_ep3);
        ys4 = repmat(wvals,1,n) .* ys3;


        end
%
