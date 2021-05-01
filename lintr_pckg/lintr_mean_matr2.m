        function [xmean,as2_mean,yback,fracs] = lintr_mean_matr2(ys,as_uniq,...
            ias,nas,m,m2,n)
%
%        Computes the least squares estimator for the mean, for general
%        rectangular transformations.
%

        icounts = zeros(1,nas);
        for i=1:nas
%
        icounts(i) = length(find(ias == i));
    end

        fracs = icounts / n;
        chk0 = sum(fracs) - 1

        as2_mean = zeros(m,m);
        yback = zeros(m,n);

        for i=1:nas
%
        as2_mean = as2_mean + fracs(i)*as_uniq(:,:,i)'*as_uniq(:,:,i);
    end

        for i=1:n
%
        yback(:,i) = as_uniq(:,:,ias(i))'*ys(:,i);
    end


        [u2,s2] = eig(as2_mean);
        s2=diag(s2);
%
        as2_inv = u2 * diag(1./s2) * u2';
        xmean = as2_inv*mean(yback,2);

%%%        chk0 = norm(as2_inv*as2_mean - eye(m))

        end
%
