        function [ys3,var_ep3] = lintr_gen2std_color(ys,as,m,n,k,var_ep)
%
%        converts general observations A*X + ep to standard spiked model
%        with colored noise by backprojecting and normalizing
%
%        returns the transformed data, and the new noise variances
%
        ys2 = zeros(m,n);
        ys3 = zeros(m,n);

        as2 = conj(as) .* as;
        as2_mean = mean(as2,2);

%
%        make back-projected data ys2, normalized backprojected ys3
%
        ys2 = conj(as).*ys;
        ys3 = ys2 ./ repmat(as2_mean,1,n);
%
        var_ep2 = as2_mean .* var_ep;
        var_ep3 = var_ep2 ./ as2_mean.^2;


        end
%
