        function [ys3,wmat,winv,as2_mean,as2_inv,spec] = lintr_gen2whi_matr(ys,as,...
            m,n,cov_ep);
%
%        backproject the data
%
        ys2 = zeros(m,n);
        for i=1:n
%
        ys2(:,i) = as(:,:,i)'*ys(:,i);
    end

%
%        construct normalization and whitening matrices
%
        as2 = zeros(m,m,n);
        covs_ep2 = zeros(m,m,n);
        for i=1:n
%
        as2(:,:,i) = as(:,:,i)'*as(:,:,i);
        covs_ep2(:,:,i) = as(:,:,i)'*cov_ep*as(:,:,i);
    end
        as2_mean = mean(as2,3);
        cov_ep2 = mean(covs_ep2,3);

        [uep2,sep2] = eig(cov_ep2);
        sep2=diag(sep2);
        wmat = uep2 * diag(1./sqrt(sep2)) * uep2';
        winv = uep2 * diag(sqrt(sep2)) * uep2';

        [ua2,sa2] = eig(as2_mean);
        spec=diag(sa2);
        as2_inv = ua2 * diag(1./spec) * ua2';

%
%        whiten the data and shrink using Donoho-Gavish
%
        ys3 = wmat * ys2;



        end
%
%
%
%
%
