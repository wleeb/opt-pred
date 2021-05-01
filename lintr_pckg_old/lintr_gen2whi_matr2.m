        function [ys3,wmat,winv,as2_mean,as2_inv,spec] = lintr_gen2whi_matr2(ys,...
            as_uniq,nas,ias,m,n,cov_ep)

        icounts = zeros(1,nas);
        for i=1:nas
%
        icounts(i) = length(find(ias == i));
    end

        fracs = icounts / n;

%
%
%        backproject the data
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys2(:,i) = as_uniq(:,:,ias(i))'*ys(:,i);
    end

%
%        construct normalization and whitening matrices
%
        cov_ep2 = zeros(m,m);
        as2_mean = zeros(m,m);

        for i=1:nas
%
        as2_mean = as2_mean + fracs(i)*as_uniq(:,:,i)'*as_uniq(:,:,i);
        cov_ep2 = cov_ep2 + fracs(i)*as_uniq(:,:,i)'*cov_ep*as_uniq(:,:,i);
    end

        [uep2,sep2] = eig(cov_ep2);
        sep2=diag(sep2);
        wmat = uep2 * diag(1./sqrt(sep2)) * uep2';
        winv = uep2 * diag(sqrt(sep2)) * uep2';

        [ua2,sa2] = eig(as2_mean);
        spec=diag(sa2);
        as2_inv = ua2 * diag(1./spec) * ua2';

        ys3 = wmat * ys2;

        end
%
%
%
%
%
