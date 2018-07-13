        function [xmean,as2_mean,yback] = lintr_mean_matr(ys,as_cube,m,m2,n)
%
%        Computes the least squares estimator for the mean, for general
%        rectangular transformations.
%
        as2_mean = zeros(m,m);
        yback = zeros(m,n);

        for i=1:n
%
        as2_mean = as2_mean + as_cube(:,:,i)'*as_cube(:,:,i) / n;
        yback(:,i) = as_cube(:,:,i)'*ys(:,i);
    end

        [u2,s2] = eig(as2_mean);
        s2=diag(s2);
%
        as2_inv = u2 * diag(1./s2) * u2';
        xmean = as2_inv*mean(yback,2);

%%%        chk0 = norm(as2_inv*as2_mean - eye(m))

        end
%
%
%
%
%
