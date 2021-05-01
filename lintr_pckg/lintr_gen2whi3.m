        function [ys4,wvals,as2_mean] = lintr_gen2whi3(ys,as,m,n,k,sig)
%
%        converts general observations A*X + ep to standard spiked model
%        with whitened noise by backprojecting, normalizing, and whitening
%
%        returns the whitened data, and the diagonal whitening transformation
%
%
%         . . . make normalization matrix
%

        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);

%
%        make back-projected data ys2 and normalized backprojected ys3
%
        ys2 = conj(as) .* ys;
        ys3 = ys2  ./ repmat(as2_mean,1,n);
%
        var_ep2 = sig^2 * as2_mean;
        var_ep3 = var_ep2 ./ as2_mean.^2;
        wvals = 1./sqrt(var_ep3);
%
%        apply transformation to whiten the noise
%
        ys4 = repmat(wvals,1,n) .* ys3;

        end
%
