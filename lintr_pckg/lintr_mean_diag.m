        function [xmean,as2_mean,yback] = lintr_mean_diag(ys,as,m,n)
%
%        Computes the least squares estimator for the mean, for diagonal
%        transformations.
%

        as2_mean = mean(abs(as).^2,2);
        yback = conj(as) .* ys;
        xmean = mean(yback,2) ./ as2_mean;


        return

%
%        the slow, direct way:
%
        as2_mean = zeros(m,1);

%
%        backproject the data, and compute normalization matrix
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys2(:,i) = conj(as(:,i)) .* ys(:,i);
        as2_mean = as2_mean + abs(as(:,i)).^2 / n;
    end

        xmean = mean(ys2,2) ./ as2_mean;

        end
%
